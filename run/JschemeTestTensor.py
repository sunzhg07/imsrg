#!/usr/bin/env python3

import pyIMSRG

emax = 2
ms = pyIMSRG.ModelSpace(emax, 'He6', 'He6')
ut = pyIMSRG.UnitTest(ms)
passed = True
#ModelSpace& modelspace, int jrank, int tz, int parity, int particle_rank, int hermitian
Jx, px, Tx, rankx, hx = 0, 0, 0, 1, -1
Jy, py, Ty, ranky, hy = 0, 0, 0, 1, +1
X = ut.RandomOp(ms, Jx, Tx, px, rankx, hx)
Y = ut.RandomOp(ms, Jy, Ty, py, ranky, hy)
X.MakeReduced();
Y.MakeReduced();
print(X.GetJRank())
passed &= ut.Mscheme_Test_comm220tt(X,Y)

print('passed?', passed)
prof = pyIMSRG.IMSRGProfiler()
prof.PrintTimes()
exit(not passed)
