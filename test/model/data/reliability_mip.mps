NAME          RELIABILITY_MIP
OBJSENSE
 MINIMIZE
ROWS
 N  OBJ
 G  XMIN
 G  YMIN
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    X         OBJ                1.0   XMIN               1.0
    Y         OBJ                1.0   YMIN               1.0
    MARKER                 'MARKER'                 'INTEND'
RHS
    RHS       XMIN               0.5   YMIN               0.5
BOUNDS
 LI BND       X                  0.0
 UI BND       X                  2.0
 LI BND       Y                  0.0
 UI BND       Y                  2.0
ENDATA
