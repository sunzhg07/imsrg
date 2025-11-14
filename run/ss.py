#!/usr/bin/env python3


from pyIMSRG import *

emax=2
ref='He4'
val='p-shell'

ms = ModelSpace(emax,ref,val)

ut = UnitTest(ms)

H = ut.RandomOp(ms, 0,0,0,2,+1 )
O = ut.RandomOp(ms, 0,0,0,2,-1 )

Res = Operator(ms,0,0,0,3)
Res.ThreeBody.SetMode('pn')
Commutator.SetIMSRG3Onlyvvv(True)
Commutator.comm223ss(H,O,Res)

print('O,H, 2b norms',O.TwoBodyNorm(),H.TwoBodyNorm())
print('Res 3b norm= ',Res.ThreeBodyNorm())
Res=Res*0.1

rw = ReadWrite()

#Res.PrintThreeBody()
rw.WriteValence3body(Res.ThreeBody, 'test.ini')


prof = IMSRGProfiler()
prof.PrintTimes()

