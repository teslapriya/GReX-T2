import sqlite3
import logging


def create_tables(con: sqlite3.Connection):
    """Create tables in the SQLite database if they don't already exist"""
    cur = con.cursor()
    cur.execute(
        """
    CREATE TABLE IF NOT EXISTS candidate (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        dm REAL NOT NULL,
        snr REAL NOT NULL,
        mjd REAL NOT NULL,
        boxcar INTEGER NOT NULL,
        sample INTEGER NOT NULL
    ) STRICT;"""
    )
    cur.execute(
        """
    CREATE TABLE IF NOT EXISTS cluster (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        centroid INTEGER NOT NULL,
        injection INTEGER,
        FOREIGN KEY (centroid) REFERENCES candidate (id),
        FOREIGN KEY (injection) REFERENCES injection (id)
    ) STRICT;"""
    )
    cur.execute(
        """
    CREATE TABLE IF NOT EXISTS cluster_member (
        candidate INTEGER PRIMARY KEY,
        cluster INTEGER NOT NULL,
        FOREIGN KEY (candidate) REFERENCES candidate (id),
        FOREIGN KEY (cluster) REFERENCES cluster (id)
    ) WITHOUT ROWID;"""
    )
    cur.execute(
        """
    CREATE TABLE IF NOT EXISTS trigger (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        cand_name TEXT,
        cluster INTEGER NOT NULL,
        FOREIGN KEY (cluster) REFERENCES cluster (id)
    );"""
    )
    con.commit()


def connect_and_create(path: str) -> sqlite3.Connection:
    """Connect to a SQLite database and create our tables if they don't exist"""
    con = sqlite3.connect(path)
    create_tables(con)
    logging.info("Successfully setup and connected to SQLite database")
    return con


def is_injection(mjd: float, con: sqlite3.Connection) -> bool:
    """Run a SQL query to see if T0 performed an injection near candidate time"""
    OFFSET = 5.78704e-5  # 5 Seconds in days
    cur = con.cursor()
    cur.execute(
        "SELECT COUNT(*) FROM injection WHERE mjd BETWEEN ? AND ?",
        mjd - OFFSET / 2,
        mjd + OFFSET / 2,
    ).fetchall()
    return cur.fetchone()[0] == 1
