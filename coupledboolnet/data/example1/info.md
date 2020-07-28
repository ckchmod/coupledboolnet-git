Hierarchical Differentiation of Myeloid Progenitors Is Encoded in the Transcription Factor Network

- set of 11 central myeloid transcription factors, known to orchestrate the respective differentiation decisions.

GATA-2 = GATA-2 AND NOT (GATA-1 AND FOG-1) AND NOT (PU.1)
GATA-1 = (GATA-1 OR GATA-2 OR Fli-1) AND NOT (PU.1)
FOG-1 = GATA-1
EKLF = GATA-1 AND NOT (Fli-1)
Fli-1 = GATA-1 AND NOT (EKLF)
SCL = GATA-1 AND NOT (PU.1)
C/EBPα = C/EBPα AND NOT (GATA-1 AND FOG-1 AND SCL)
PU.1 = (C/EBPα OR PU.1) AND NOT (GATA-1 OR GATA-2)
cJun = PU.1 AND NOT (Gfi-1)
EgrNab = (PU.1 AND cJun) AND NOT (Gfi-1)
Gfi-1 = C/EBPα AND NOT (EgrNab)


x_1 = x_1 AND NOT (x_2 AND x_3) AND NOT (x_8)
x_2 = (x_2 OR x_1 OR x_5) AND NOT (x_8)
x_3 = x_2
x_4 = x_2 AND NOT (x_5)
x_5 = x_2 AND NOT (x_4)
x_6 = x_2 AND NOT (x_8)
x_7 = x_7 AND NOT (x_2 AND x_3 AND x_6)
x_8 = (x_7 OR x_8) AND NOT (x_2 OR x_1)
x_9 = x_8 AND NOT (x_11)
x_10 = (x_8 AND x_9) AND NOT (x_11)
x_11 = x_7 AND NOT (x_10)

     f_1 | f_2 | f_3 | f_4 | f_5 | f_6 | f_7 | f_8 | f_9  | f_10 | f_11   
_______________________________________________________________________
1  | 0   | 0   | 0   | 0   | 0   | 0   |  0  |  0  |  0   |  0   |  0
2  | 0   | 0   | 1   | 0   | 0   | 0   |  0  |  0  |  0   |  0   |  0
3  | 0   | 1   |     | 1   | 1   | 1   |  0  |  0  |  1   |  0   |  1
4  | 0   | 0   |     | 0   | 0   | 0   |  0  |  0  |  0   |  0   |  0
5  | 0   | 1                           |  0  |  1         |  0
6  | 0   | 0                           |  0  |  0         |  0
7  | 0   | 1                           |  0  |  0         |  1
8  | 0   | 0                           |  0  |  0         |  0
9  | 1   | 1                           |  1  |  1
10 | 0   | 0                           |  1  |  0
11 | 1   | 1                           |  1  |  0
12 | 0   | 0                           |  1  |  0
13 | 1   | 1                           |  1  |  1
14 | 0   | 0                           |  1  |  0
15 | 0   | 1                           |  1  |  0
16 | 0   | 0                           |  0  |  0
-----------------------------------------------------------------------
     x_1 | x_2 | x_2 | x_2 | x_2 | x_2 | x_7 | x_7 | x_8  | x_8  | x_7   
     x_2 | x_1 |     | x_5 | x_4 | x_8 | x_2 | x_8 | x_11 | x_9  | x_10
     x_3 | x_5 |                       | x_3 | x_2 |      | x_11 |
     x_8 | x_8 |                       | x_6 | x_1 | 