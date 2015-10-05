__author__ = 'ruben'
__doc__ = 'Query the sqlite database of Buszaki hc database to extract cells, cellType, region. It outputs a txt file'

import sqlite3 as sqlite

tablesToIgnore = ["sqlite_sequence"]
outputFilename = None


def Print(msg):
    if (outputFilename != None):
        outputFile = open(outputFilename, 'a')
        print >> outputFile, msg
        outputFile.close()
    else:
        print msg


def Describe(dbFile):
    connection = sqlite.connect(dbFile)
    cursor = connection.cursor()

    Print("TableName\tColumns\tRows\tCells")

    totalTables = 0
    totalColumns = 0
    totalRows = 0
    totalCells = 0

    # Get List of Tables:
    tableListQuery = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY Name"
    cursor.execute(tableListQuery)
    tables = map(lambda t: t[0], cursor.fetchall())

    for table in tables:

        if (table in tablesToIgnore):
            continue

        columnsQuery = "PRAGMA table_info(%s)" % table
        cursor.execute(columnsQuery)
        numberOfColumns = len(cursor.fetchall())

        rowsQuery = "SELECT Count() FROM %s" % table
        cursor.execute(rowsQuery)
        numberOfRows = cursor.fetchone()[0]

        numberOfCells = numberOfColumns * numberOfRows

        Print("%s\t%d\t%d\t%d" % (table, numberOfColumns, numberOfRows, numberOfCells))

        totalTables += 1
        totalColumns += numberOfColumns
        totalRows += numberOfRows
        totalCells += numberOfCells

    Print("")
    Print("Number of Tables:\t%d" % totalTables)
    Print("Total Number of Columns:\t%d" % totalColumns)
    Print("Total Number of Rows:\t%d" % totalRows)
    Print("Total Number of Cells:\t%d" % totalCells)

    return cursor, connection


if __name__ == "__main__":
    dbFile = '/media/bigdata/hc-3/hc3-metadata-tables/hc3-tables.db'
    basep = dbFile.split('/hc3-tables.db')
    cursor, connection = Describe(dbFile)

    linear_query = "SELECT cellType, region FROM cell WHERE topdir='ec013.15'"
    cursor.execute(linear_query)

    cell_type = cursor.fetchall()

    print '{} Cells found'.format(len(cell_type))

    with open(basep[0] + "/isIntern.txt", "w") as text_f:
        text_f.write("cellId\tisIntern?\tregion\n")
        for idx, i in enumerate(cell_type):
            t = 1 if i[0] == 'i' else 0
            text_f.write("{}\t{}\t{}\n".format(idx + 1, t, i[1]))

    print('File saved')
    cursor.close()
    connection.close()
