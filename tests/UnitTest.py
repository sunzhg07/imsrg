import pyIMSRG 

ms=pyIMSRG.ModelSpace(2,'He4','He4')
ut=pyIMSRG.UnitTest(ms)
passed=ut.SanityCheck()
op=pyIMSRG.Operator(ms)
hole=next(iter(ms.holes))
op.SetOneBody(hole, hole, 1.0)
do_no=op.DoNormalOrdering()
undo_no=op.UndoNormalOrdering()
try:
    assert passed
    assert abs(do_no.ZeroBody - undo_no.ZeroBody) > 1.0e-12
    print("All tests passed")
except AssertionError:  
    print("Tests failed")
    raise AssertionError
