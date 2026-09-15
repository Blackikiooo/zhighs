NAME          PRESOLVE_MIP
ROWS
 N  OBJ
 G  DEMAND
 E  FIXROW
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    X         OBJ                1.0   DEMAND             1.0
    Y         OBJ                1.0   DEMAND             1.0
    Z         OBJ                5.0   DEMAND             1.0
    Z         FIXROW             1.0
    MARKER                 'MARKER'                 'INTEND'
RHS
    RHS       DEMAND             2.0
    RHS       FIXROW             1.0
BOUNDS
 BV BND       X
 BV BND       Y
 UP BND       Z                  1.0
 LO BND       Z                  1.0
ENDATA
