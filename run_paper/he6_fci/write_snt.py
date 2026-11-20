import sys
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "..", "build")); sys.path.insert(0, os.path.join(HERE, ".."))
from pyIMSRG import ModelSpace, Operator, ReadWrite
from he8_io import hs_path
ms = ModelSpace(3, "He4", "p-shell"); ms.SetHbarOmega(16)
H = Operator(ms, 0, 0, 0, 3); ReadWrite().ReadOperator(H, hs_path("He6", 16, 3, val="p-shell"))
print("He4-core ZeroBody", H.ZeroBody)
Hvac = H.UndoNormalOrdering()
msfci = ModelSpace(3, "vacuum", "FCI"); msfci.SetHbarOmega(16)
Hvac.SetModelSpace(msfci)
print("vacuum ZeroBody", Hvac.ZeroBody)
ReadWrite().WriteTokyo(Hvac, sys.argv[1], "")
