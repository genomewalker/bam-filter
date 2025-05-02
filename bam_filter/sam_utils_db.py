import duckdb
import logging
import os

log = logging.getLogger("my_logger")


def load_sam_file(con, file_path):
    """Load entire SAM file as raw lines into DuckDB table"""
    # Drop existing tables if they exist to avoid conflicts
    con.execute("DROP TABLE IF EXISTS raw_lines")
    con.execute("DROP TABLE IF EXISTS header")
    con.execute("DROP TABLE IF EXISTS alignments")
    con.execute("DROP TABLE IF EXISTS qname_index")
    con.execute("DROP TABLE IF EXISTS rname_index")

    log.info("Reading SAM file header lines...")
    # First statement: Load the raw data - process in smaller chunks
    raw_query = """
    CREATE TABLE raw_lines AS
    SELECT 
        row_number() OVER () as line_id,
        line 
    FROM read_csv_auto(
        ?, 
        delim='\n',           -- Read line by line
        header=FALSE, 
        names=['line'],       -- Single column for entire line
        sample_size=100000,   -- Sample size for schema detection
        all_varchar=TRUE,     -- Ensure all lines are read as strings
        parallel=TRUE,        -- Enable parallel reading
        buffer_size=262144    -- Use larger buffer for better I/O
    );
    """
    # Execute query - no configuration changes after this point
    con.execute(raw_query, [file_path])

    log.info("Extracting header lines...")
    header_query = """
    CREATE TABLE header AS 
    SELECT line as header_line
    FROM raw_lines
    WHERE line LIKE '@%';
    """
    con.execute(header_query)

    log.info("Processing alignment data - optimized for performance...")
    # Optimized alignment processing - split line once and use array indexing
    alignments_query = """
    CREATE TABLE alignments AS 
    WITH split_data AS (
        SELECT 
            line, 
            string_split(line, '\t') AS parts
        FROM raw_lines
        WHERE NOT line LIKE '@%'
    )
    SELECT
        parts[1] AS qname,
        TRY_CAST(parts[2] AS INTEGER) AS flag,
        NULLIF(parts[3], '*') AS rname,
        TRY_CAST(parts[4] AS INTEGER) AS pos,
        TRY_CAST(parts[5] AS INTEGER) AS mapq,
        NULLIF(parts[6], '*') AS cigar,
        NULLIF(parts[7], '*') AS rnext,
        TRY_CAST(parts[8] AS INTEGER) AS pnext,
        TRY_CAST(parts[9] AS INTEGER) AS tlen,
        NULLIF(parts[10], '*') AS seq,
        NULLIF(parts[11], '*') AS qual,
        -- Extract tags more efficiently
        CASE 
            WHEN len(parts) > 11 THEN list_slice(parts, 12, len(parts))
            ELSE list_value(NULL) 
        END AS tags
    FROM split_data;
    """

    # Execute query directly without trying to change any database settings
    con.execute(alignments_query)

    # Get counts before dropping raw_lines
    raw_count = con.execute("SELECT COUNT(*) FROM raw_lines").fetchone()[0]
    header_count = con.execute("SELECT COUNT(*) FROM header").fetchone()[0]
    alignment_count = con.execute("SELECT COUNT(*) FROM alignments").fetchone()[0]

    # Create summary/index tables for faster queries
    log.info("Creating optimized index tables for faster queries...")

    # Create qname summary table - useful for read-based operations
    con.execute(
        """
    CREATE TABLE qname_index AS
    SELECT 
        qname,
        COUNT(*) as alignment_count,
        array_agg(DISTINCT rname) as references,
        MIN(mapq) as min_mapq,
        MAX(mapq) as max_mapq
    FROM alignments
    GROUP BY qname;
    
    -- Create index on qname for faster lookups
    CREATE INDEX idx_qname ON qname_index(qname);
    """
    )

    # Create rname summary table - useful for reference-based operations
    con.execute(
        """
    CREATE TABLE rname_index AS
    SELECT 
        rname,
        COUNT(*) as read_count,
        COUNT(DISTINCT qname) as unique_reads,
        AVG(mapq) as avg_mapq
    FROM alignments
    WHERE rname IS NOT NULL
    GROUP BY rname;
    
    -- Create index on rname for faster lookups
    CREATE INDEX idx_rname ON rname_index(rname);
    """
    )

    # Track number of index entries
    qname_index_count = con.execute("SELECT COUNT(*) FROM qname_index").fetchone()[0]
    rname_index_count = con.execute("SELECT COUNT(*) FROM rname_index").fetchone()[0]
    log.info(
        f"Created index tables: {qname_index_count} qnames, {rname_index_count} rnames"
    )

    # Drop raw data immediately to free memory
    con.execute("DROP TABLE raw_lines")

    # Log results
    log.info(
        f"Processed {raw_count} total lines: {header_count} header lines and {alignment_count} alignment records"
    )
    return f"Header and alignment data successfully loaded ({alignment_count} records)"


# Add new utility functions for working with the indexes


def get_reads_by_reference(con, reference_name, min_mapq=0):
    """Get all reads mapping to a specific reference"""
    query = """
    SELECT a.*
    FROM alignments a
    WHERE a.rname = ? AND a.mapq >= ?
    """
    return con.execute(query, [reference_name, min_mapq])


def get_references_by_read(con, read_name):
    """Get all references a specific read maps to"""
    query = """
    SELECT r.*, a.mapq, a.flag
    FROM alignments a
    JOIN rname_index r ON a.rname = r.rname
    WHERE a.qname = ?
    """
    return con.execute(query, [read_name])


def get_multimapped_reads(con, min_alignments=2):
    """Get reads that map to multiple references"""
    query = """
    SELECT q.qname, q.alignment_count, q.references
    FROM qname_index q
    WHERE array_length(q.references) >= ?
    """
    return con.execute(query, [min_alignments])


def format_sam_output(con):
    """Format alignments for SAM output"""
    query = """
    SELECT 
        qname,
        CAST(flag AS VARCHAR),
        COALESCE(rname, '*'),
        CAST(pos AS VARCHAR),
        CAST(mapq AS VARCHAR),
        COALESCE(cigar, '*'),
        COALESCE(rnext, '*'),
        CAST(pnext AS VARCHAR),
        CAST(tlen AS VARCHAR),
        COALESCE(seq, '*'),
        COALESCE(qual, '*'),
        ARRAY_TO_STRING(tags, '\t')
    FROM alignments
    """
    return query
