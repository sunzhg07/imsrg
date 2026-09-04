#include "IMSRG.hh"
#include "FactorizedDoubleCommutator_eths.hh"
#include "version.hh"
#include <Python.h>
#include <sstream>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//  Orbit MS_GetOrbit(ModelSpace& self, int i){ return self.GetOrbit(i);};
//  size_t MS_GetOrbitIndex_Str(ModelSpace& self, std::string s){ return
//  self.GetOrbitIndex(s);}; TwoBodyChannel MS_GetTwoBodyChannel(ModelSpace&
//  self, int ch){return self.GetTwoBodyChannel(ch);};

//  double TB_GetTBME_J(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int
//  c, int d){return self.GetTBME_J(j_bra,j_ket,a,b,c,d);}; double
//  TB_GetTBME_J_norm(TwoBodyME& self,int j_bra, int j_ket, int a, int b, int c,
//  int d){return self.GetTBME_J_norm(j_bra,j_ket,a,b,c,d);};

//  size_t TBCGetLocalIndex(TwoBodyChannel& self, int p, int q){ return
//  self.GetLocalIndex( p, q);};

//  void ArmaMatPrint( arma::mat& self){ self.print();};
//  void OpSetOneBodyME( Operator& self, int i, int j, double
//  v){self.OneBody(i,j) = v;};

//  void MS_SetRef(ModelSpace& self, std::string str){ self.SetReference(
//  str);}; void MS_SetRef(ModelSpace& self, const std::set<index_t>& ref){
//  self.SetReference( ref);};

//  Operator HF_GetNormalOrderedH(HartreeFock& self){ return
//  self.GetNormalOrderedH();}; Operator HF_GetNormalOrderedH(HartreeFock& self,
//  int particle_rank=2){ return self.GetNormalOrderedH(particle_rank);};

// BOOST_PYTHON_MODULE(pyIMSRG)
// PYBIND11_PLUGIN(pyIMSRG)
PYBIND11_MODULE(pyIMSRG, m) {
  m.doc() = "python bindings for IMSRG code";

  py::class_<Orbit>(m, "Orbit")
      .def(py::init<>())
      .def_readwrite("n", &Orbit::n)
      .def_readwrite("l", &Orbit::l)
      .def_readwrite("j2", &Orbit::j2)
      .def_readwrite("tz2", &Orbit::tz2)
      .def_readwrite("occ", &Orbit::occ)
      .def_readwrite("cvq", &Orbit::cvq)
      .def_readwrite("index", &Orbit::index);

  py::class_<TwoBodyChannel>(m, "TwoBodyChannel")
      .def(py::init<>())
      .def("GetNumberKets", &TwoBodyChannel::GetNumberKets)
      //      .def("GetLocalIndex",&TBCGetLocalIndex)
      .def("GetLocalIndex", [](TwoBodyChannel &self, int p,
                               int q) { return self.GetLocalIndex(p, q); })
      .def("GetKetIndex", &TwoBodyChannel::GetKetIndex)
      .def("GetKet", [](TwoBodyChannel &self, int i) { return self.GetKet(i); })
      .def("GetKetIndex_pp",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_pp();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_hh",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_hh();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_ph",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_ph();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_cc",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_cc();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_vc",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_vc();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_qc",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_qc();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_vv",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_vv();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_qv",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_qv();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def("GetKetIndex_qq",
           [](TwoBodyChannel &self) {
             auto &x = self.GetKetIndex_qq();
             std::vector<size_t> v(x.begin(), x.end());
             return v;
           })
      .def_readwrite("J", &TwoBodyChannel::J)
      .def_readwrite("parity", &TwoBodyChannel::parity)
      .def_readwrite("Tz", &TwoBodyChannel::Tz);

      py::class_<TwoBodyChannel_CC>(m, "TwoBodyChannel_CC")
          .def(py::init<>())
          .def("GetNumberKets", &TwoBodyChannel_CC::GetNumberKets)
          .def("GetKet", [](TwoBodyChannel_CC &self, int i)
               { return self.GetKet(i); })
          .def("GetLocalIndex",
               [](TwoBodyChannel_CC &self, int p, int q)
               { return self.GetLocalIndex(p, q); })
          .def_readwrite("J", &TwoBodyChannel_CC::J)
          .def_readwrite("parity", &TwoBodyChannel_CC::parity)
          .def_readwrite("Tz", &TwoBodyChannel_CC::Tz);

      py::class_<ThreeBodyChannel>(m, "ThreeBodyChannel")
          .def(py::init<>())
          .def("GetNumber3bKets", &ThreeBodyChannel::GetNumber3bKets)
          .def("GetLocalIndex", &ThreeBodyChannel::GetLocalIndex, py::arg("p"), py::arg("q"), py::arg("r"), py::arg("Jpq"))
          .def("GetKet", [](ThreeBodyChannel &self, int i)
               { return self.GetKet(i); })
          .def_readwrite("twoJ", &ThreeBodyChannel::twoJ)
          .def_readwrite("parity", &ThreeBodyChannel::parity)
          .def_readwrite("twoTz", &ThreeBodyChannel::twoTz);

  py::class_<Ket>(m, "Ket")
      .def(py::init<Orbit &, Orbit &>())
      .def_readwrite("p", &Ket::p)
      .def_readwrite("q", &Ket::q);
  py::class_<Ket3>(m, "Ket3")
      .def(py::init<Orbit &, Orbit &, Orbit &>())
      .def_readwrite("p", &Ket3::p)
      .def_readwrite("q", &Ket3::q)
      .def_readwrite("r", &Ket3::r)
      .def_readwrite("Jpq", &Ket3::Jpq);

      py::class_<ModelSpace>(m, "ModelSpace")
          .def(py::init<>())
          .def(py::init<const ModelSpace &>())
          .def(py::init<int, const std::string &>(), py::arg("emax"), py::arg("reference"))
          .def(py::init<int, const std::string &, const std::string &>(), py::arg("emax"), py::arg("reference"), py::arg("valence"))
          .def(py::init<int, std::vector<std::string>, std::vector<std::string>>(), py::arg("emax"), py::arg("hole_list"), py::arg("valence_list"))
          .def(py::init<int, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>(), py::arg("emax"), py::arg("hole_list"), py::arg("core_list"), py::arg("valence_list"))
          .def("SetHbarOmega", &ModelSpace::SetHbarOmega)
          .def("SetTargetMass", &ModelSpace::SetTargetMass)
          .def("SetTargetZ", &ModelSpace::SetTargetZ)
          .def(
              "AddOrbit", [](ModelSpace &self, int n, int l, int j2, int tz2, double occ, int cvq)
              { self.AddOrbit(n, l, j2, tz2, occ, cvq); },
              py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"), py::arg("occ"), py::arg("cvq"))
          .def("SetupKets", &ModelSpace::SetupKets)
          .def("Setup3bKets", &ModelSpace::Setup3bKets)
          .def("SetOcc", &ModelSpace::SetOcc, py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"), py::arg("occ"))
          .def("SetOccNAT", &ModelSpace::SetOccNAT, py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"), py::arg("occ_nat"))
          .def("SetEmax", &ModelSpace::SetEmax)
          .def("SetE2max", &ModelSpace::SetE2max)
          .def("SetE3max", &ModelSpace::SetE3max)
          .def("SetdE3max", &ModelSpace::SetdE3max)
          .def("SetLmax", &ModelSpace::SetLmax)
          .def("SetEmaxUnocc", &ModelSpace::SetEmaxUnocc)
          .def("SetEmax3Body", &ModelSpace::SetEmax3Body)
          .def("FindEFermi", &ModelSpace::FindEFermi)
          .def("GetHbarOmega", &ModelSpace::GetHbarOmega)
          .def("GetTargetMass", &ModelSpace::GetTargetMass)
          .def("GetTargetZ", &ModelSpace::GetTargetZ)
          .def("GetAref", &ModelSpace::GetAref)
               .def("GetAcore", &ModelSpace::GetAcore)
          .def("GetZref", &ModelSpace::GetZref)
          .def("GetNumberOrbits", &ModelSpace::GetNumberOrbits)
          .def("GetNumberKets", &ModelSpace::GetNumberKets)
          .def("GetNumberTwoBodyChannels", &ModelSpace::GetNumberTwoBodyChannels)
          .def("GetNumberTwoBodyChannels_CC", &ModelSpace::GetNumberTwoBodyChannels_CC)
          .def("GetNumberThreeBodyChannels", &ModelSpace::GetNumberThreeBodyChannels)
          .def("GetEmax",  &ModelSpace::GetEmax)
          .def("GetE2max", &ModelSpace::GetE2max)
          .def("GetE3max", &ModelSpace::GetE3max)
          //      .def("GetOrbit", &MS_GetOrbit)
          .def("GetOrbit", [](ModelSpace &self, int i)
               { return self.GetOrbit(i); })
          .def("GetKet", [](ModelSpace &self, int i)
               { return self.GetKet(i); })
          .def("GetTwoBodyChannelIndex", &ModelSpace::GetTwoBodyChannelIndex)
          .def("GetTwoBodyChannel", [](ModelSpace &self, int ch)
               { return self.GetTwoBodyChannel(ch); })
          .def("GetTwoBodyChannel_CC", [](ModelSpace &self, int ch)
               { return self.GetTwoBodyChannel_CC(ch); })
          .def("GetThreeBodyChannel", &ModelSpace::GetThreeBodyChannel)
          .def("GetThreeBodyChannelIndex", &ModelSpace::GetThreeBodyChannelIndex, py::arg("twoJ"), py::arg("parity"), py::arg("twoTz"))
          .def("Index2String", &ModelSpace::Index2String)
          .def("ResetFirstPass", &ModelSpace::ResetFirstPass)
          //      .def("SetReference", &MS_SetRef)
          .def("SetReference", [](ModelSpace &self, const std::set<index_t> &ref)
               { self.SetReference(ref); })
          .def("SetReferenceStr", [](ModelSpace &self, std::string s)
               { self.SetReference(s); })
          .def("SetReferenceOcc", [](ModelSpace &self, std::map<index_t,double> &ref)
               { self.SetReference(ref); })
          .def("Init_occ_from_file", &ModelSpace::Init_occ_from_file)
          .def("InitSingleSpecies", &ModelSpace::InitSingleSpecies)
          .def(
              "GetOrbitIndex", [](ModelSpace &self, int n, int l, int j, int tz)
              { return self.GetOrbitIndex(n, l, j, tz); },
              py::arg("n"), py::arg("l"), py::arg("j2"), py::arg("tz2"))
          .def(
              "GetOrbitIndex_fromString", [](ModelSpace &self, std::string s)
              { return self.GetOrbitIndex(s); },
              py::arg("orbstring"))
          .def(
              "GetOneBodyChannels", [](ModelSpace &self, int l, int j, int tz)
              { return self.OneBodyChannels.at({l, j, tz}); },
              py::arg("l"), py::arg("j2"), py::arg("tz2"))
          //      .def("GetOrbitIndex_fromString", &MS_GetOrbitIndex_Str)
          .def("PreCalculateSixJ", &ModelSpace::PreCalculateSixJ)
          .def("PreCalculateNineJ", &ModelSpace::PreCalculateNineJ)
          .def("PreCalculateMoshinsky",&ModelSpace::PreCalculateMoshinsky)
          .def("GetMoshinsky",&ModelSpace::GetMoshinsky)
          .def("GetSixJ",&ModelSpace::GetSixJ)
          .def("GetNineJ",&ModelSpace::GetNineJ)
          .def("NineJHash",&ModelSpace::NineJHash)
//          .def("NineJUnHash",&ModelSpace::NineJUnHash)
          .def("NineJUnHash",[](ModelSpace &self, uint64_t key){ uint64_t k1,k2,k3,k4,k5,k6,k7,k8,k9; self.NineJUnHash(key,k1,k2,k3,k4,k5,k6,k7,k8,k9); return py::make_tuple(k1,k2,k3,k4,k5,k6,k7,k8,k9);  }     )
          .def("SetScalarFirstPass", &ModelSpace::SetScalarFirstPass)
          .def("SetScalar3bFirstPass", &ModelSpace::SetScalar3bFirstPass)
          .def("ClearVectors", &ModelSpace::ClearVectors)
          .def("Print", &ModelSpace::Print)
//          .def("Print_CC", &ModelSpace::Print_CC)
          .def_readwrite("holes", &ModelSpace::holes)
          .def_readwrite("particles", &ModelSpace::particles)
          .def_readwrite("core", &ModelSpace::core)
          .def_readwrite("valence", &ModelSpace::valence)
          .def_readwrite("qspace", &ModelSpace::qspace)
          .def_readwrite("all_orbits", &ModelSpace::all_orbits);

      py::class_<Operator>(m, "Operator")
          .def(py::init<>())
          .def(py::init<ModelSpace &>())
          .def(py::init<Operator &>())
          .def(py::init<ModelSpace &, int, int, int, int>(), py::arg("modelspace"), py::arg("j_rank"), py::arg("t_rank"), py::arg("parity"), py::arg("particle_rank"))
          .def(py::self += py::self)
          .def(py::self + py::self)
          .def(py::self -= Operator())
          .def(py::self - Operator())
          .def(-py::self)
          .def(py::self *= double())
          .def(py::self * double())
          .def(double() * py::self)
          .def(py::self /= double())
          .def(py::self / double())
          .def(py::self += double())
          .def(py::self + double())
          .def(py::self -= double())
          .def(py::self - double())
          .def_readwrite("ZeroBody", &Operator::ZeroBody)
          .def_readwrite("OneBody", &Operator::OneBody)
          .def_readwrite("TwoBody", &Operator::TwoBody)
          .def_readwrite("ThreeBody", &Operator::ThreeBody)
          .def("GetOneBody", &Operator::GetOneBody, py::arg("i"), py::arg("j"))
          .def("SetOneBody", &Operator::SetOneBody, py::arg("i"), py::arg("j"), py::arg("MatEl"))
          .def("GetTwoBody", &Operator::GetTwoBody, py::arg("ch_bra"), py::arg("ch_ket"), py::arg("ibra"), py::arg("iket"))
          .def("SetTwoBody", &Operator::SetTwoBody)
          .def("GetTwoBodyDimension", &Operator::GetTwoBodyDimension)
          .def("ScaleOneBody", &Operator::ScaleOneBody)
          .def("ScaleTwoBody", &Operator::ScaleTwoBody)
          .def("EraseOneBody", &Operator::EraseOneBody)
          .def("EraseTwoBody", &Operator::EraseTwoBody)
          .def("EraseThreeBody", &Operator::EraseThreeBody)
          .def("DoNormalOrdering", &Operator::DoNormalOrdering)
          .def("DoNormalOrderingCore", &Operator::DoNormalOrderingCore)
          .def("DoNormalOrderingFilledValence", &Operator::DoNormalOrderingFilledValence)
          .def("UndoNormalOrdering", &Operator::UndoNormalOrdering)
          .def("UndoNormalOrderingCore", &Operator::UndoNormalOrderingCore)
          .def("ReNormalOrderCore", &Operator::ReNormalOrderCore)
          .def("DoNormalOrdering", &Operator::UndoNormalOrdering)
          .def("SetModelSpace", &Operator::SetModelSpace)
          .def("GetModelSpace", &Operator::GetModelSpace,
               py::return_value_policy::reference)
          .def("Truncate", &Operator::Truncate)
          .def("DoIsospinAveraging", &Operator::DoIsospinAveraging)
          .def("Norm", &Operator::Norm)
          .def("OneBodyNorm", &Operator::OneBodyNorm)
          .def("OneLegNorm", &Operator::OneLegNorm)
          .def("TwoBodyNorm", &Operator::TwoBodyNorm)
          .def("ThreeBodyNorm", &Operator::ThreeBodyNorm)
          .def("SetHermitian", &Operator::SetHermitian)
          .def("SetAntiHermitian", &Operator::SetAntiHermitian)
          .def("SetNonHermitian", &Operator::SetNonHermitian)
          .def("IsHermitian", &Operator::IsHermitian)
          .def("IsAntiHermitian", &Operator::IsAntiHermitian)
          .def("IsReduced", &Operator::IsReduced)
          .def("PrintOneBody", &Operator::PrintOneBody)
          .def("PrintTwoBody", [](Operator &self)
               { self.PrintTwoBody(); })
          .def("PrintTwoBody_ch", [](Operator &self, int ch)
               { self.PrintTwoBody(ch); })
          .def("PrintTwoBody_chch", [](Operator &self, int ch_bra, int ch_ket)
               { self.PrintTwoBody(ch_bra, ch_ket); })
          .def("PrintThreeBody", &Operator::PrintThreeBody )
          //      .def("PrintTwoBody_ch", &Operator::PrintTwoBody)
          .def("MakeReduced", &Operator::MakeReduced)
          .def("MakeNotReduced", &Operator::MakeNotReduced)
          .def("MakeNormalized", &Operator::MakeNormalized)
          .def("MakeUnNormalized", &Operator::MakeUnNormalized)
          .def("GetParticleRank", &Operator::GetParticleRank)
          .def("SetParticleRank", &Operator::SetParticleRank)
          .def("GetJRank", &Operator::GetJRank)
          .def("GetTRank", &Operator::GetTRank)
          .def("GetParity", &Operator::GetParity)
          .def("GetNumberLegs", &Operator::GetNumberLegs)
          .def("SetNumberLegs", &Operator::SetNumberLegs, py::arg("legs"))
          .def("SetQSpaceOrbit", &Operator::SetQSpaceOrbit, py::arg("q"))
          .def("GetQSpaceOrbit", &Operator::GetQSpaceOrbit)
          .def("IsNumberConserving", &Operator::IsNumberConserving)
          .def("ThreeLegNorm", &Operator::ThreeLegNorm)
          .def("EraseThreeLeg", &Operator::EraseThreeLeg)
          .def_readwrite("ThreeLeg", &Operator::ThreeLeg)
//          .def("GetE3max", &Operator::GetE3max)
//          .def("SetE3max", &Operator::SetE3max)
          .def("PrintTimes", &Operator::PrintTimes)
          .def("Size", &Operator::Size)
          .def("MakeNormalized", &Operator::MakeNormalized)
          .def("MakeUnNormalized", &Operator::MakeUnNormalized)
          .def("GetOneBodyChannel", &Operator::GetOneBodyChannel, py::arg("l"), py::arg("j2"), py::arg("tz2"))
          //      .def("SetOneBodyME", &OpSetOneBodyME)
          .def("SetOneBodyME", [](Operator &self, int i, int j, double v)
               { self.OneBody(i, j) = v; })
          .def("GetMP2_Energy", &Operator::GetMP2_Energy)
          .def("GetMP2_3BEnergy", &Operator::GetMP2_Energy)
          .def("GetMP3_Energy", &Operator::GetMP3_Energy)
          .def("GetPPHH_Ladders", &Operator::GetPPHH_Ladders)
          .def(
              "ReadBinary", [](Operator &self, std::string fname)
              { std::ifstream ifs(fname,std::ios::binary);  self.ReadBinary(ifs); },
              py::arg("filename"))
          .def(
              "WriteBinary", [](Operator &self, std::string fname)
              { std::ofstream ofs(fname,std::ios::binary);  self.WriteBinary(ofs); },
              py::arg("filename"))
          //      .def("IsospinProject", &Operator::IsospinProject)
          ;

      py::class_<ThreeLegME>(m, "ThreeLegME")
          .def(py::init<>())
          .def("GetME", &ThreeLegME::GetME, py::arg("ch"), py::arg("a"), py::arg("b"), py::arg("c"))
          .def("GetME_norm", &ThreeLegME::GetME_norm, py::arg("ch"), py::arg("a"), py::arg("b"), py::arg("c"))
          .def("GetME_J", &ThreeLegME::GetME_J, py::arg("J"), py::arg("a"), py::arg("b"), py::arg("c"))
          .def("SetME", &ThreeLegME::SetME, py::arg("ch"), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("me"))
          .def("AddToME", &ThreeLegME::AddToME, py::arg("ch"), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("me"))
          .def("AddToME_J", &ThreeLegME::AddToME_J, py::arg("J"), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("me"))
          .def("Norm", &ThreeLegME::Norm)
          .def("Erase", &ThreeLegME::Erase)
          .def("Allocate", &ThreeLegME::Allocate)
          .def("GetMatrix", [](ThreeLegME &self, size_t ch) -> arma::mat &
               { return self.GetMatrix(ch); },
               py::arg("ch"), py::return_value_policy::reference_internal);

  py::class_<arma::mat>(m, "ArmaMat")
      .def(py::init<>())
      .def(
          "zeros",
          [](arma::mat &self, int nrows, int ncols) {
            self.zeros(nrows, ncols);
          },
          py::arg("nrows"), py::arg("ncols"))
      .def("Print", [](arma::mat &self) { self.print(); }) //   &ArmaMatPrint)
      .def("__str__",
           [](arma::mat &self) {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           }) //   &ArmaMatPrint)
      .def(
          "save", [](arma::mat &self, std::string fname) { self.save(fname); },
          py::arg("filename"))
      .def(
          "load", [](arma::mat &self, std::string fname) { self.load(fname); },
          py::arg("filename"))
      //      .def("t", &arma::mat::t) // transpose
      .def("t",
           [](arma::mat &self) {
             arma::mat x = self.t();
             return x;
           }) // transpose
      .def(py::self *= double())
      //      .def(py::self * double())
      //      .def(double() * py::self)
      //      .def(double() * py::self, [](double x, arma::mat& self){arma::mat
      //      out = x * self; return out;} )
      .def(py::self /= double())
      .def(py::self / double())
      //      .def(py::self += ArmaMat())
      //      .def(py::self + ArmaMat())
      //      .def(py::self -= ArmaMat())
      //      .def(py::self - ArmaMat())
      .def(
          "__mul__",
          [](const arma::mat &A, const arma::mat &B) {
            arma::mat C = A * B;
            return C;
          },
          py::is_operator())
      .def(
          "__mul__",
          [](const arma::mat &B, float A) {
            arma::mat C = A * B;
            return C;
          },
          py::is_operator())
      //      .def("__mul__", [](float A, const arma::mat& B){arma::mat C = A *
      //      B; return C;}, py::is_operator() )
      .def(
          "__add__",
          [](const arma::mat &A, const arma::mat &B) {
            arma::mat C = A + B;
            return C;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const arma::mat &A, const arma::mat &B) {
            arma::mat C = A - B;
            return C;
          },
          py::is_operator())
      .def(
          "__call__",
          [](arma::mat &self, const int i, const int j) { return &self(i, j); },
          py::is_operator())
      .def(
          "Set",
          [](arma::mat &self, const int i, const int j, double x) {
            self(i, j) = x;
          },
          py::arg("i"), py::arg("j"), py::arg("matel"))
      .def("Getn_rows", [](arma::mat &self) { return self.n_rows; })
      .def("Getn_cols", [](arma::mat &self) { return self.n_cols; })
      .def("Schur_Prod",
           [](arma::mat &self, arma::mat &other) {
             arma::mat out = self % other;
             return out;
           })
      .def("Norm", [](arma::mat &self) { return arma::norm(self, "fro"); })
      .def("trace",
           [](arma::mat &self) {
             double t = arma::trace(self);
             return t;
           })
      .def("sum", [](arma::mat &self) {
        double s = arma::accu(self);
        return s;
      });

  py::class_<TwoBodyME>(m, "TwoBodyME")
      .def(py::init<>())
      //      .def("GetTBME_J", TB_GetTBME_J)
      //      .def("GetTBME_J_norm", TB_GetTBME_J_norm)
      .def("GetTBME_J",
           [](TwoBodyME &self, int Jbra, int Jket, int a, int b, int c, int d) {
             return self.GetTBME_J(Jbra, Jket, a, b, c, d);
           })
      .def("GetTBME_J_norm",
           [](TwoBodyME &self, int Jbra, int Jket, int a, int b, int c, int d) {
             return self.GetTBME_J_norm(Jbra, Jket, a, b, c, d);
           })
      .def(
          "GetTBMEmonopole",
          [](TwoBodyME &self, int a, int b, int c, int d) {
            return self.GetTBMEmonopole(a, b, c, d);
          },
          py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"))
      .def("GetTBME_norm",
           [](TwoBodyME &self, int ch_bra, int ch_ket, int a, int b, int c,
              int d) { return self.GetTBME_norm(ch_bra, ch_ket, a, b, c, d); })
      .def(
          "GetTBMEmonopole_norm",
          [](TwoBodyME &self, int a, int b, int c, int d) {
            return self.GetTBMEmonopole_norm(a, b, c, d);
          },
          py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"))
      .def(
          "GetChannelMatrix",
          [](TwoBodyME &self, int J, int p, int Tz) {
            size_t ch = self.modelspace->GetTwoBodyChannelIndex(J, p, Tz);
            return self.GetMatrix(ch, ch);
          },
          py::arg("J"), py::arg("parity"), py::arg("Tz"))
      .def("PrintAll",
           [](TwoBodyME &self) {
             for (auto &it : self.MatEl) {
               if (it.second.n_rows > 0) {
                 std::cout << it.first[0] << " " << it.first[1] << std::endl
                           << it.second << std::endl;
               };
             };
           })
      .def("PrintMatrix", &TwoBodyME::PrintMatrix, py::arg("ch_bra"),
           py::arg("ch_ket"))
      .def("Erase", &TwoBodyME::Erase)
      .def("GetTBMEnorm_chij",
           [](TwoBodyME &self, int ch_bra, int ch_ket, size_t ibra,
              size_t iket) {
             return self.GetTBME_norm(ch_bra, ch_ket, ibra, iket);
           })
      .def("SetTBME_chij",
           [](TwoBodyME &self, int ch_bra, int ch_ket, size_t ibra, size_t iket,
              double tbme) { self.SetTBME(ch_bra, ch_ket, ibra, iket, tbme); })
      .def("Norm", &TwoBodyME::Norm)
      .def(py::self *= double())
      .def(double() * py::self)
      .def(py::self * double())
      .def(py::self + TwoBodyME())
      .def(py::self += TwoBodyME())
      .def(py::self - TwoBodyME())
      .def(py::self -= TwoBodyME());

  //   py::class_<ThreeBodyME>(m,"ThreeBodyME")
  //      .def(py::init<>())
  //      .def("SetME", &ThreeBodyME::SetME)
  //      .def("GetME", &ThreeBodyME::GetME)
  //      .def("GetME_pn", &ThreeBodyME::GetME_pn)
  //      .def("RecouplingCoefficient",&ThreeBodyME::RecouplingCoefficient)
  //      .def_readonly_static("ABC",&ThreeBodyME::ABC)
  //      .def_readonly_static("BCA",&ThreeBodyME::BCA)
  //      .def_readonly_static("CAB",&ThreeBodyME::CAB)
  //      .def_readonly_static("ACB",&ThreeBodyME::ACB)
  //      .def_readonly_static("CBA",&ThreeBodyME::CBA)
  //      .def_readonly_static("BAC",&ThreeBodyME::BAC)
  //   ;

  //   py::class_<ThreeBodyMEpn>(m,"ThreeBodyMEpn")
  py::class_<ThreeBodyME>(m, "ThreeBodyME")
      .def(py::init<>())
      //      .def("SetME", &ThreeBodyMEpn::SetME)
      //      .def("GetME", &ThreeBodyME::GetME)
      .def(
          "GetME_iso",
          [](ThreeBodyME &self, int Jab, int Jde, int twoJ, int tab, int tde,
             int twoTabc, int twoTdef, int a, int b, int c, int d, int e,
             int f) {
            return self.GetME_iso(Jab, Jde, twoJ, tab, tde, twoTabc, twoTdef, a,
                                  b, c, d, e, f);
          },
          py::arg("Jab"), py::arg("Jde"), py::arg("twoJ"), py::arg("tab"),
          py::arg("tde"), py::arg("twoTabc"), py::arg("twoTdef"), py::arg("a"),
          py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"), py::arg("f"))
      .def("SetME_pn", &ThreeBodyME::SetME_pn)
      // .def("GetME_pn", &ThreeBodyME::GetME_pn)
      .def(
          "GetME_pn",
          [](ThreeBodyME &self, int Jab_in, int Jde_in, int twoJ, int a, int b,
             int c, int d, int e, int f) {
            return self.GetME_pn(Jab_in, Jde_in, twoJ, a, b, c, d, e, f);
          },
          py::arg("Jab_in"), py::arg("Jde_in"), py::arg("twoJ"), py::arg("a"),
          py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"), py::arg("f"))
      .def(
          "GetME_pn_tensor",
          [](ThreeBodyME &self, int Jab_in, int j0, int Jde_in, int j1, int a,
             int b, int c, int d, int e, int f) {
            return self.GetME_pn(Jab_in, j0, Jde_in, j1, a, b, c, d, e, f);
          },
          py::arg("Jab_in"), py::arg("j0"), py::arg("Jde_in"), py::arg("j1"),
          py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"),
          py::arg("f"))
      .def("SetME_pn_ch", &ThreeBodyME::SetME_pn_ch)
      .def("GetME_pn_ch", &ThreeBodyME::GetME_pn_ch)
      .def("GetME_pn_no2b", &ThreeBodyME::GetME_pn_no2b)
      .def("Get_ch_start_keys", [](const ThreeBodyME &self) {
          // Return list of (ch_bra, ch_ket) pairs that are actually stored
          std::vector<std::pair<size_t,size_t>> keys;
          for (auto &kv : self.Get_ch_start())
              keys.emplace_back(kv.first.ch_bra, kv.first.ch_ket);
          return keys;
      })
      .def("RecouplingCoefficient", &ThreeBodyME::RecouplingCoefficient)
      .def("TransformToPN", &ThreeBodyME::TransformToPN)
      .def("SwitchToPN_and_discard", &ThreeBodyME::SwitchToPN_and_discard)
      //      .def("Print",&ThreeBodyME::Print)
      //      .def("PrintAll",&ThreeBodyME::PrintAll)
      .def("Erase", &ThreeBodyME::Erase)
      .def("SetMode", &ThreeBodyME::SetMode)
      .def("IsAllocated", &ThreeBodyME::IsAllocated)
      .def("IsHermitian", &ThreeBodyME::IsHermitian)
      .def("Is_PN_Mode", &ThreeBodyME::Is_PN_Mode)
      .def("ReadFile", &ThreeBodyME::ReadFile, py::arg("string_inputs"),
           py::arg("int_inputs"))
      .def(py::self += ThreeBodyME(), py::is_operator())
      .def(py::self *= double())
      //      .def_readonly_static("ABC",&ThreeBodyME::ABC)
      //      .def_readonly_static("BCA",&ThreeBodyME::BCA)
      //      .def_readonly_static("CAB",&ThreeBodyME::CAB)
      //      .def_readonly_static("ACB",&ThreeBodyME::ACB)
      //      .def_readonly_static("CBA",&ThreeBodyME::CBA)
      //      .def_readonly_static("BAC",&ThreeBodyME::BAC)
      ;

  py::class_<ReadWrite>(m, "ReadWrite")
      .def(py::init<>())
      .def("ReadTBME_Oslo", &ReadWrite::ReadTBME_Oslo)
      .def("ReadTBME_OakRidge", &ReadWrite::ReadTBME_OakRidge,
           py::arg("spname"), py::arg("tbmename"), py::arg("H"),
           py::arg("tbme_format") = "ascii")
      .def("ReadBareTBME_Jason", &ReadWrite::ReadBareTBME_Jason)
      .def("ReadBareTBME_Navratil", &ReadWrite::ReadBareTBME_Navratil)
      .def("ReadBareTBME_Darmstadt", &ReadWrite::ReadBareTBME_Darmstadt,
           py::arg("filename"), py::arg("H"), py::arg("e1max"),
           py::arg("e2max"), py::arg("lmax"))
      .def("Read_Darmstadt_3body", &ReadWrite::Read_Darmstadt_3body,
           py::arg("filename"), py::arg("H"), py::arg("e1max"),
           py::arg("e2max"), py::arg("e3max"))
      .def("ReadOperator2b_Miyagi", &ReadWrite::ReadOperator2b_Miyagi,
           py::arg("filename"), py::arg("ms"))
#ifndef NO_HDF5
      .def("Read3bodyHDF5", &ReadWrite::Read3bodyHDF5)
#endif
      .def("Write_me2j", &ReadWrite::Write_me2j)
      .def("Write_me2j_gz", &ReadWrite::Write_me2j_gz)
      .def("Write_me3j", &ReadWrite::Write_me3j)
      .def("WriteTBME_Navratil", &ReadWrite::WriteTBME_Navratil)
      .def("WriteNuShellX_sps", &ReadWrite::WriteNuShellX_sps, py::arg("op"),
           py::arg("filename"))
      .def("WriteNuShellX_int", &ReadWrite::WriteNuShellX_int, py::arg("op"),
           py::arg("filename"))
      .def("WriteNuShellX_op", &ReadWrite::WriteNuShellX_op, py::arg("op"),
           py::arg("filename"))
      .def("ReadNuShellX_int", &ReadWrite::ReadNuShellX_int, py::arg("op"),
           py::arg("filename"))
      .def("ReadNuShellX_int_iso", &ReadWrite::ReadNuShellX_int_iso,
           py::arg("op"), py::arg("filename"))
      .def("WriteAntoine_int", &ReadWrite::WriteAntoine_int)
      .def("WriteAntoine_input", &ReadWrite::WriteAntoine_input)
      .def("WriteOperator", &ReadWrite::WriteOperator)
      .def("WriteOperatorHuman", &ReadWrite::WriteOperatorHuman)
      .def("ReadOperator", &ReadWrite::ReadOperator)
      .def("ReadOperatorHuman", &ReadWrite::ReadOperatorHuman)
      .def("CompareOperators", &ReadWrite::CompareOperators)
      .def("WriteOneBody_Simple", &ReadWrite::WriteOneBody_Simple)
      .def("ReadOneBody_Takayuki", &ReadWrite::ReadOneBody_Takayuki)
      .def("ReadTwoBody_Takayuki", &ReadWrite::ReadTwoBody_Takayuki)
      .def("WriteOneBody_Takayuki", &ReadWrite::WriteOneBody_Takayuki)
      .def("WriteTwoBody_Takayuki", &ReadWrite::WriteTwoBody_Takayuki)
      .def("WriteTensorOneBody", &ReadWrite::WriteTensorOneBody)
      .def("WriteTensorTwoBody", &ReadWrite::WriteTensorTwoBody)
      .def("WriteTokyo", &ReadWrite::WriteTokyo, py::arg("op"),
           py::arg("filename"), py::arg("mode"))
      .def("WriteTensorTokyo", &ReadWrite::WriteTensorTokyo,
           py::arg("filename"), py::arg("op"))
      .def(
          "ReadTokyo",
          [](ReadWrite &self, std::string s, Operator &op) {
            self.ReadTokyo(s, op);
          },
          py::arg("file_in"), py::arg("op"))
      .def(
          "ReadTensorTokyo",
          [](ReadWrite &self, std::string s, Operator &op) {
            self.ReadTensorTokyo(s, op);
          },
          py::arg("file_in"), py::arg("op"))
      .def("WriteOneBody_Oslo", &ReadWrite::WriteOneBody_Oslo)
      .def("WriteTwoBody_Oslo", &ReadWrite::WriteTwoBody_Oslo)
      .def("SetCoMCorr", &ReadWrite::SetCoMCorr)
      .def("ReadTwoBodyEngel", &ReadWrite::ReadTwoBodyEngel)
      .def("ReadOperator_Nathan", &ReadWrite::ReadOperator_Nathan)
      .def("ReadTensorOperator_Nathan", &ReadWrite::ReadTensorOperator_Nathan)
      .def("ReadRelCMOpFromJavier", &ReadWrite::ReadRelCMOpFromJavier)
      .def("Set3NFormat", &ReadWrite::Set3NFormat)
      .def("WriteDaggerOperator", &ReadWrite::WriteDaggerOperator)
      .def("ReadJacobi3NFiles", &ReadWrite::ReadJacobi3NFiles)
      .def("WriteValence3body", &ReadWrite::WriteValence3body)
      .def("ReadValence3body",  &ReadWrite::ReadValence3body)
      .def("SetScratchDir", &ReadWrite::SetScratchDir)
      .def("GetScratchDir", &ReadWrite::GetScratchDir)
      .def("CopyFile", &ReadWrite::CopyFile, py::arg("filein"),
           py::arg("fileout"))
      .def("ReadDarmstadt_2bodyRel", &ReadWrite::ReadDarmstadt_2bodyRel)
      .def("ReadH2_2body", &ReadWrite::ReadH2_2body)
      .def("Read2bCurrent_Navratil", &ReadWrite::Read2bCurrent_Navratil,
           py::arg("filename"), py::arg("Op"))
      //      .def("WriteOmega",&ReadWrite::WriteOmega,
      //      py::arg("basename"),py::arg("scratch_dir"),py::arg("nOmegas"))
      ;

  py::class_<HartreeFock>(m, "HartreeFock")
      .def(py::init<Operator &>())
      .def("Solve", &HartreeFock::Solve)
      .def("TransformToHFBasis", &HartreeFock::TransformToHFBasis)
      .def("GetHbare", &HartreeFock::GetHbare)
      //      .def("GetNormalOrderedH",&HF_GetNormalOrderedH)
      //      .def("GetNormalOrderedH",&HF_GetNormalOrderedH,
      //      py::arg("particle_rank")=2 )
      .def(
          "GetNormalOrderedH",
          [](HartreeFock &self, int pRank) {
            return self.GetNormalOrderedH(pRank);
          },
          py::arg("particle_rank") = 2)
      .def(
          "GetNormalOrderedH_Cin",
          [](HartreeFock &self, arma::mat &C, int pRank) {
            return self.GetNormalOrderedH(C, pRank);
          },
          py::arg("C"), py::arg("particle_rank") = 2)
      .def("GetOmega", &HartreeFock::GetOmega)
      .def("PrintSPE", &HartreeFock::PrintSPE)
      .def("PrintSPEandWF", &HartreeFock::PrintSPEandWF)
      .def("GetRadialWF_r", &HartreeFock::GetRadialWF_r)
      .def("GetHFPotential", &HartreeFock::GetHFPotential)
      .def("GetAverageHFPotential", &HartreeFock::GetAverageHFPotential)
      .def("GetValence3B", &HartreeFock::GetValence3B)
      .def("FreeVmon", &HartreeFock::FreeVmon)
      .def("UpdateDensityMatrix", &HartreeFock::UpdateDensityMatrix)
      .def("UpdateF", &HartreeFock::UpdateF)
      .def("BuildMonopoleV", &HartreeFock::BuildMonopoleV)
      .def("CalcEHF", &HartreeFock::CalcEHF)
      .def("PrintEHF", &HartreeFock::PrintEHF)
      .def("FillLowestOrbits", &HartreeFock::FillLowestOrbits)
      .def("DiscardNO2Bfrom3N", &HartreeFock::DiscardNO2Bfrom3N)
      .def("FreezeOccupations", &HartreeFock::FreezeOccupations)
      .def("UnFreezeOccupations", &HartreeFock::UnFreezeOccupations)
      .def_static("Vmon3Hash", &HartreeFock::Vmon3Hash)
      // Modifying arguments which were passed by reference causes trouble in
      // python, so instead we bind a lambda function and return a tuple
      .def_static("Vmon3UnHash",
                  [](uint64_t key) {
                    int a, b, c, d, e, f;
                    HartreeFock::Vmon3UnHash(key, a, b, c, d, e, f);
                    return std::make_tuple(a, b, c, d, e, f);
                  })
      .def_readonly("EHF", &HartreeFock::EHF)
      .def_readonly("F", &HartreeFock::F)     // Fock matrix
      .def_readonly("rho", &HartreeFock::rho) // density matrix
                                              //      .def_readonly("C",&HartreeFock::C)
                                              //      // Unitary transformation
      .def_readwrite("C", &HartreeFock::C) // Unitary transformation
      .def_readwrite("Vmon3_keys", &HartreeFock::Vmon3_keys)
      .def_readwrite("Vmon3", &HartreeFock::Vmon3);

  py::class_<HFMBPT, HartreeFock>(m, "HFMBPT")
      .def(py::init<Operator &>())
      .def("UseNATOccupations", &HFMBPT::UseNATOccupations)
      .def("GetNaturalOrbitals", &HFMBPT::GetNaturalOrbitals)
      .def("TransformHOToNATBasis", &HFMBPT::TransformHOToNATBasis)
      .def("TransformHFToNATBasis", &HFMBPT::TransformHFToNATBasis)
      .def("GetNormalOrderedHNAT", &HFMBPT::GetNormalOrderedHNAT)
      .def("PrintSPEandWF", &HFMBPT::PrintSPEandWF)
      .def_readwrite("C_HO2NAT", &HFMBPT::C_HO2NAT) // Unitary transformation
      .def_readwrite("C_HF2NAT", &HFMBPT::C_HF2NAT) // Unitary transformation
      ;

  // Define which overloaded version of IMSRGSolver::Transform I want to expose
  //   Operator (IMSRGSolver::*Transform_ref)(Operator&) =
  //   &IMSRGSolver::Transform;

  py::class_<IMSRGSolver>(m, "IMSRGSolver")
      .def(py::init<Operator &>())
      .def("Solve", &IMSRGSolver::Solve)
      //      .def("Transform",Transform_ref)
      .def("Transform",
           [](IMSRGSolver &self, Operator &op) { return self.Transform(op); })
      .def("InverseTransform", &IMSRGSolver::InverseTransform)
      .def("SetFlowFile", &IMSRGSolver::SetFlowFile)
      .def("SetMethod", &IMSRGSolver::SetMethod)
      .def("SetEtaCriterion", &IMSRGSolver::SetEtaCriterion)
      .def("SetDs", &IMSRGSolver::SetDs)
      .def("SetdOmega", &IMSRGSolver::SetdOmega)
      .def("SetOmegaNormMax", &IMSRGSolver::SetOmegaNormMax)
      .def("SetSmax", &IMSRGSolver::SetSmax)
      .def("SetDsmax", &IMSRGSolver::SetDsmax)
      .def("SetHin", &IMSRGSolver::SetHin)
      .def("SetODETolerance", &IMSRGSolver::SetODETolerance)
      .def("Reset", &IMSRGSolver::Reset)
      .def("SetGenerator", &IMSRGSolver::SetGenerator)
      .def("SetOnly2bEta",
           [](IMSRGSolver &self, bool tf) {
             self.GetGenerator().SetOnly2bEta(tf);
           })
      .def("SetDenominatorCutoff", &IMSRGSolver::SetDenominatorCutoff)
      .def("SetDenominatorDelta", &IMSRGSolver::SetDenominatorDelta)
      .def("SetDenominatorDeltaOrbit", &IMSRGSolver::SetDenominatorDeltaOrbit)
      .def("SetDenominatorPartitioning",
           &IMSRGSolver::SetDenominatorPartitioning) // Can be Epstein_Nesbet
                                                     // (default) or
                                                     // Moller_Plesset
      .def("GetSystemDimension", &IMSRGSolver::GetSystemDimension)
      // .def("GetOmega", &IMSRGSolver::GetOmega)
      .def("GetOmega", py::overload_cast<int>(&IMSRGSolver::GetOmega),
           py::arg("index"), "Get an Operator at a specific index")
      .def("GetOmega", py::overload_cast<>(&IMSRGSolver::GetOmega),
           "Get the entire deque of Operators")
      .def("SetOmega", &IMSRGSolver::SetOmega, py::arg("index"),
           py::arg("Omega"))
      //      .def("GetH_s",&IMSRGSolver::GetH_s,return_value_policy<reference_existing_object>())
      .def("GetH_s", &IMSRGSolver::GetH_s)
      .def("SetH_s", &IMSRGSolver::SetH_s)
      .def("GetS", &IMSRGSolver::GetS)
      .def("SetMagnusAdaptive", &IMSRGSolver::SetMagnusAdaptive)
      .def("SetReadWrite", &IMSRGSolver::SetReadWrite)
      .def("SetHunterGatherer", &IMSRGSolver::SetHunterGatherer)
      //          .def("SetPerturbativeTriples",
      //          &IMSRGSolver::SetPerturbativeTriples)
      //          .def("GetPerturbativeTriples",
      //          &IMSRGSolver::GetPerturbativeTriples)
      //          .def("CalculatePerturbativeTriples",
      //          &IMSRGSolver::CalculatePerturbativeTriples)
      .def("CalculatePerturbativeTriples",
           py::overload_cast<>(&IMSRGSolver::CalculatePerturbativeTriples))
      .def("CalculatePerturbativeTriples",
           py::overload_cast<Operator &>(
               &IMSRGSolver::CalculatePerturbativeTriples))
      .def("AddOperator", &IMSRGSolver::AddOperator)
      .def("GetOperator", &IMSRGSolver::GetOperator)
      .def("EstimateBCHError", &IMSRGSolver::EstimateBCHError)
      .def("UpdateEta", &IMSRGSolver::UpdateEta)
      .def("GetNOmegaWritten", &IMSRGSolver::GetNOmegaWritten)
      .def("GetOmegaSize", &IMSRGSolver::GetOmegaSize)
      //      .def("GetScratchDir",[](IMSRGSolver& self){ return
      //      self.rw->GetScratchDir();} )
      .def("GetScratchDir", [](IMSRGSolver &self) { return self.scratchdir; })
      .def("FlushOmegaToScratch", &IMSRGSolver::FlushOmegaToScratch)
      .def_readwrite("generator", &IMSRGSolver::generator)
      .def_readwrite("Eta", &IMSRGSolver::Eta)
      .def_readwrite("n_omega_written",
                     &IMSRGSolver::n_omega_written) // I'm not sure I like just
                                                    // directly exposing this...
      .def("SetOnly1bEta", [](IMSRGSolver &self, bool tf) {
        self.GetGenerator().SetOnly1bEta(tf);
      });

  py::class_<IMSRGSolverPV, IMSRGSolver>(m, "IMSRGSolverPV")
      .def(py::init<Operator &, Operator &>())
      .def_readwrite("Etapv", &IMSRGSolverPV::Etapv)
      .def("Solve_RK4", &IMSRGSolverPV::Solve_flow_RK4_PV)
      .def("Solve_magnus_euler", &IMSRGSolverPV::Solve_magnus_euler_PV)
      .def("AddOperatorPV", &IMSRGSolverPV::AddOperatorPV)
      .def("GetOperatorPV", &IMSRGSolverPV::GetOperatorPV)
      .def("GetVPT_s", &IMSRGSolverPV::GetVPT_s)
      .def("SetGeneratorPV", &IMSRGSolverPV::SetGeneratorPV)
      .def("SetOnly1bEta",
           [](IMSRGSolverPV &self, bool tf) {
             self.GetGeneratorPV().SetOnly1bEta(tf);
           })
      .def("Transform", [](IMSRGSolverPV &self, Operator &op, Operator &opPV) {
        return self.Transform(op, opPV);
      });

  py::class_<Generator>(m, "Generator")
      .def(py::init<>())
      .def("SetType", &Generator::SetType, py::arg("gen_type"))
      .def("SetDenominatorPartitioning", &Generator::SetDenominatorPartitioning,
           py::arg("Moller_Plessett or Epstein_Nesbet"))
      .def("SetUseIsospinAveraging", &Generator::SetUseIsospinAveraging,
           py::arg("tf"))
      .def("Update", &Generator::Update, py::arg("H"), py::arg("Eta"))
      .def("GetHod_SingleRef", &Generator::GetHod_SingleRef, py::arg("H"))
      .def("GetHod", &Generator::GetHod, py::arg("H"));

  py::class_<GeneratorPV, Generator>(m, "GeneratorPV")
      .def(py::init<>())
      .def("SetType", &Generator::SetType, py::arg("gen_type"))
      .def("Update", &GeneratorPV::Update, py::arg("H"), py::arg("V"),
           py::arg("Eta"), py::arg("Etapv"));

  py::class_<IMSRGProfiler>(m, "IMSRGProfiler")
      .def(py::init<>())
      .def("PrintTimes", &IMSRGProfiler::PrintTimes)
      .def("PrintCounters", &IMSRGProfiler::PrintCounters)
      .def("PrintAll", &IMSRGProfiler::PrintAll)
      .def("PrintMemory", &IMSRGProfiler::PrintMemory)
      .def("Clear", &IMSRGProfiler::Clear);

  py::class_<Jacobi3BME>(m, "Jacobi3BME")
      .def(py::init<>())
      .def(py::init<int, int, int, int, int>())
      .def("GetDimensionAS", &Jacobi3BME::GetDimensionAS)
      .def("GetDimensionNAS", &Jacobi3BME::GetDimensionNAS)
      .def("GetMatElAS", &Jacobi3BME::GetMatElAS)
      .def("GetMatElNAS", &Jacobi3BME::GetMatElNAS)
      .def("SetEmax", &Jacobi3BME::SetEmax)
      .def("SetE2max", &Jacobi3BME::SetE2max)
      .def("SetE3max", &Jacobi3BME::SetE3max)
      .def("ComputeNAS_MatrixElements", &Jacobi3BME::ComputeNAS_MatrixElements)
      .def("GetLabMatEl", &Jacobi3BME::GetLabMatEl)
      .def("TestReadTcoeffNavratil", &Jacobi3BME::TestReadTcoeffNavratil)
      //          .def("GetV3mon_all", &Jacobi3BME::GetV3mon_all)
      ;

  py::module Commutator = m.def_submodule("Commutator", "Commutator namespace");
  Commutator.def("Commutator", &Commutator::Commutator);
  Commutator.def("CommutatorScalarScalar", &Commutator::CommutatorScalarScalar);
  Commutator.def("CommutatorScalarTensor", &Commutator::CommutatorScalarTensor);
  Commutator.def("CommutatorScalarDagger", &Commutator::CommutatorScalarDagger);
  Commutator.def("comm211sd", &Commutator::comm211sd);
  Commutator.def("comm231sd", &Commutator::comm231sd);
  Commutator.def("comm413_233sd", &Commutator::comm413_233sd);
  Commutator.def("comm433_pp_hh_431sd", &Commutator::comm433_pp_hh_431sd);
  Commutator.def("comm433sd_ph", &Commutator::comm433sd_ph);
  Commutator.def("SetUseIMSRG3", &Commutator::SetUseIMSRG3);
  Commutator.def("SetUseIMSRG3N7", &Commutator::SetUseIMSRG3N7);
  Commutator.def("SetUseIMSRG3N7_Tensor", &Commutator::SetUseIMSRG3N7_Tensor);
  Commutator.def("SetUseIMSRG3_Tensor", &Commutator::SetUseIMSRG3_Tensor);
  Commutator.def("TurnOnTerm", &Commutator::TurnOnTerm);
  Commutator.def("TurnOffTerm", &Commutator::TurnOffTerm);
  Commutator.def("SetThreebodyThreshold", &Commutator::SetThreebodyThreshold);
  Commutator.def("SetVerbose", &Commutator::SetVerbose, py::arg("tf"));
  Commutator.def("SetSingleThread", &Commutator::SetSingleThread,
                 py::arg("tf"));
  Commutator.def("PrintSettings", &Commutator::PrintSettings);

  // IMSRG(2) commutators
  Commutator.def("comm110ss", &Commutator::comm110ss);
  Commutator.def("comm220ss", &Commutator::comm220ss);
  Commutator.def("comm111ss", &Commutator::comm111ss);
  Commutator.def("comm121ss", &Commutator::comm121ss);
  Commutator.def("comm221ss", &Commutator::comm221ss);
  Commutator.def("comm122ss", &Commutator::comm122ss);
  Commutator.def("comm222_pp_hh_221ss", &Commutator::comm222_pp_hh_221ss);
  Commutator.def("comm222_pp_hhss", &Commutator::comm222_pp_hhss);
  Commutator.def("comm222_phss", &Commutator::comm222_phss);
  Commutator.def("comm222_phss_slower", &Commutator::comm222_phss_slower);
  // IMSRG(3) commutators
  Commutator.def("comm330ss", &Commutator::comm330ss);
  Commutator.def("comm331ss", &Commutator::comm331ss);
  Commutator.def("comm231ss", &Commutator::comm231ss);
  Commutator.def("comm132ss", &Commutator::comm132ss);
  Commutator.def("comm232ss", &Commutator::comm232ss);
  Commutator.def("comm332_ppph_hhhpss", &Commutator::comm332_ppph_hhhpss);
  Commutator.def("comm332_pphhss", &Commutator::comm332_pphhss);
  Commutator.def("comm332ss", [](Operator &X, Operator &Y, Operator &Z) {
    Commutator::comm332_ppph_hhhpss(X, Y, Z);
    Commutator::comm332_pphhss(X, Y, Z);
  });
  Commutator.def("comm223ss", &Commutator::comm223ss);
  Commutator.def("comm133ss", &Commutator::comm133ss);
  Commutator.def("comm233_pp_hhss", &Commutator::comm233_pp_hhss);
  Commutator.def("comm233_phss", &Commutator::comm233_phss);
  Commutator.def("comm333_ppp_hhhss", &Commutator::comm333_ppp_hhhss);
  Commutator.def("comm333_pph_hhpss", &Commutator::comm333_pph_hhpss);
  // scalar-tensor commutators
  Commutator.def("comm111st", &Commutator::comm111st);
  Commutator.def("comm121st", &Commutator::comm121st);
  Commutator.def("comm221st", &Commutator::comm221st);
  Commutator.def("comm122st", &Commutator::comm122st);
  Commutator.def("comm222_pp_hh_221st", &Commutator::comm222_pp_hh_221st);
  Commutator.def("comm222_phst", &Commutator::comm222_phst);
  Commutator.def("SetIMSRG3Noqqq", &Commutator::SetIMSRG3Noqqq);
  Commutator.def("SetIMSRG3Onlyvvv", &Commutator::SetIMSRG3Onlyvvv);
  Commutator.def("SetIMSRG3valence2b", &Commutator::SetIMSRG3valence2b);
  Commutator.def("Discard0bFrom3b", &Commutator::Discard0bFrom3b);
  Commutator.def("Discard1bFrom3b", &Commutator::Discard1bFrom3b);
  Commutator.def("Discard2bFrom3b", &Commutator::Discard2bFrom3b);
  Commutator.def("comm331st", &Commutator::comm331st);
  Commutator.def("comm223st", &Commutator::comm223st);
  Commutator.def("comm223tts", &Commutator::comm223tts);
  Commutator.def("comm231st", &Commutator::comm231st);
     Commutator.def("comm231tts", &Commutator::comm231tts);
  Commutator.def("comm110tts", &Commutator::comm110tts);
  Commutator.def("comm220tts", &Commutator::comm220tts);
  Commutator.def("comm111tts", &Commutator::comm111tts);
  Commutator.def("comm121tts", &Commutator::comm121tts);
  Commutator.def("comm122tts", &Commutator::comm122tts);
  Commutator.def("comm221tts", &Commutator::comm221tts);
  Commutator.def("comm222_pp_hhtts", &Commutator::comm222_pp_hhtts);
  Commutator.def("comm222_phtts", &Commutator::comm222_phtts);
  Commutator.def("comm232st", &Commutator::comm232st);
  Commutator.def("comm133st", &Commutator::comm133st);
  Commutator.def("comm132st", &Commutator::comm132st);
     Commutator.def("comm132tts", &Commutator::comm132tts);
     Commutator.def("comm232tts", &Commutator::comm232tts);

  //      Commutator.def("comm223_231_Factorization",
  //      &Commutator::comm223_231_Factorization);
  //      Commutator.def("comm223_232_Factorization",
  //      &Commutator::comm223_232_Factorization);

  //      Commutator.def("comm223_231_Factorization_slow",
  //      &Commutator::comm223_231_Factorization_slow);
  //      Commutator.def("comm223_232_Factorization_slow",
  //      &Commutator::comm223_232_Factorization_slow);

  //       BCH.def("EstimateBCHError", &BCH::EstimateBCHError); // This doesn't
  //       really work

  py::module FactorizedDoubleCommutator = Commutator.def_submodule(
      "FactorizedDoubleCommutator", "FactorizedDoubleCommutator namespace");
  FactorizedDoubleCommutator.def(
      "comm223_231", &Commutator::FactorizedDoubleCommutator::comm223_231);
  FactorizedDoubleCommutator.def(
      "comm223_232", &Commutator::FactorizedDoubleCommutator::comm223_232);
  FactorizedDoubleCommutator.def(
      "comm223_132", &Commutator::FactorizedDoubleCommutator::comm223_132);
  FactorizedDoubleCommutator.def(
       "comm223_132_cross",
       &Commutator::FactorizedDoubleCommutator::comm223_132_cross);

  FactorizedDoubleCommutator.def(
      "comm223_231_chi2b",
      &Commutator::FactorizedDoubleCommutator::comm223_231_chi2b);
  FactorizedDoubleCommutator.def(
      "comm223_231_chi1b",
      &Commutator::FactorizedDoubleCommutator::comm223_231_chi1b);
  FactorizedDoubleCommutator.def(
      "comm223_232_chi2b",
      &Commutator::FactorizedDoubleCommutator::comm223_232_chi2b);
  FactorizedDoubleCommutator.def(
      "comm223_232_chi1b",
      &Commutator::FactorizedDoubleCommutator::comm223_232_chi1b);

  //        FactorizedDoubleCommutator.def("comm223_231_slow",
  //        &Commutator::FactorizedDoubleCommutator::comm223_231_slow);
  //        FactorizedDoubleCommutator.def("comm223_232_slow",
  //        &Commutator::FactorizedDoubleCommutator::comm223_232_slow);
  //        FactorizedDoubleCommutator.def("UseSlowVersion",
  //        &Commutator::FactorizedDoubleCommutator::UseSlowVersion);
  FactorizedDoubleCommutator.def(
      "SetUse_GooseTank_1b",
      &Commutator::FactorizedDoubleCommutator::SetUse_GooseTank_1b);
  FactorizedDoubleCommutator.def(
      "SetUse_GooseTank_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_GooseTank_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_1b_Intermediates",
      &Commutator::FactorizedDoubleCommutator::SetUse_1b_Intermediates);
  FactorizedDoubleCommutator.def(
      "SetUse_2b_Intermediates",
      &Commutator::FactorizedDoubleCommutator::SetUse_2b_Intermediates);
  FactorizedDoubleCommutator.def(
      "SetUse_GooseTank_only_1b",
      &Commutator::FactorizedDoubleCommutator::SetUse_GooseTank_only_1b);
  FactorizedDoubleCommutator.def(
      "SetUse_GooseTank_only_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_GooseTank_only_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeII_1b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeII_1b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeIII_1b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeIII_1b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeI_1b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeI_1b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeIIIa_1b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeIIIa_1b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeII_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeII_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeIII_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeIII_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGI_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGI_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGII_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGII_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGIIIa_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGIIIa_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGIIIb_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGIIIb_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGIIIc_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGIIIc_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGIVa_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGIVa_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGIVb_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGIVb_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_TypeGIVc_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_TypeGIVc_2b);

  py::module FactorizedDoubleCommutator_eths = Commutator.def_submodule(
       "FactorizedDoubleCommutator_eths",
       "Tensor-focused FactorizedDoubleCommutator_eths namespace");
  FactorizedDoubleCommutator_eths.def(
       "comm223_231_st",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_231_st);
  FactorizedDoubleCommutator_eths.def(
       "comm223_231_chi1b_tensor",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_231_chi1b_tensor);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_1b_Intermediates",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_1b_Intermediates);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_2b_Intermediates",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_2b_Intermediates);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeI_1b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeI_1b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeII_1b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeII_1b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeIII_1b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeIII_1b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeIIIa_1b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeIIIa_1b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeIIIa_slow",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeIIIa_slow);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGI_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGI_2b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGII_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGII_2b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGIIIa_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGIIIa_2b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGIIIb_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGIIIb_2b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGIIIc_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGIIIc_2b);
  FactorizedDoubleCommutator_eths.def(
       "DebugChiPandyaHermiticity",
       &Commutator::FactorizedDoubleCommutator_eths::DebugChiPandyaHermiticity);
  FactorizedDoubleCommutator_eths.def(
       "DebugTensorPandyaRoundTrip",
       &Commutator::FactorizedDoubleCommutator_eths::DebugTensorPandyaRoundTrip);
  FactorizedDoubleCommutator_eths.def(
       "ForceScalarMakeNotReduced",
       &Commutator::FactorizedDoubleCommutator_eths::ForceScalarMakeNotReduced);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232_GIIIa",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232_GIIIa);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232_GIIIb",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232_GIIIb);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232_GIIIc",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232_GIIIc);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGIVa_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGIVa_2b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGIVb_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGIVb_2b);
  FactorizedDoubleCommutator_eths.def(
       "SetUse_TypeGIVc_2b",
       &Commutator::FactorizedDoubleCommutator_eths::SetUse_TypeGIVc_2b);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232_GIVa",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232_GIVa);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232_GIVb",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232_GIVb);
  FactorizedDoubleCommutator_eths.def(
       "comm223_232_GIVc",
       &Commutator::FactorizedDoubleCommutator_eths::comm223_232_GIVc);
  FactorizedDoubleCommutator.def(
      "SetUse_GT_TypeI_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_GT_TypeI_2b);
  FactorizedDoubleCommutator.def(
      "SetUse_GT_TypeIV_2b",
      &Commutator::FactorizedDoubleCommutator::SetUse_GT_TypeIV_2b);

  py::module BCH = m.def_submodule("BCH", "BCH namespace");
  BCH.def("BCH_Transform", &BCH::BCH_Transform);
  BCH.def("BCH_Product", &BCH::BCH_Product);
  BCH.def("SetUseFactorizedCorrection", &BCH::SetUseFactorizedCorrection);
  BCH.def("SetUseFactorizedCorrectionBCH_product",
          &BCH::SetUseFactorizedCorrectionBCH_product);
  BCH.def("SetUseFactorized_Correct_ZBTerm",
          &BCH::SetUseFactorized_Correct_ZBTerm);
  BCH.def("SetOnly2bOmega", &BCH::SetOnly2bOmega);
  BCH.def("SetComm223_231", &BCH::SetComm223_231);
  BCH.def("SetComm223_232", &BCH::SetComm223_232);
  BCH.def("Set_BCH_Transform_Threshold", &BCH::Set_BCH_Transform_Threshold);
  BCH.def("Set_BCH_Product_Threshold", &BCH::Set_BCH_Product_Threshold);
  BCH.def("SetBCHSkipiEq1", &BCH::SetBCHSkipiEq1);

  py::module ReferenceImplementations = m.def_submodule(
      "ReferenceImplementations", "ReferenceImplementations namespace");
  ReferenceImplementations.def("comm110ss",
                               &ReferenceImplementations::comm110ss);
  ReferenceImplementations.def("comm220ss",
                               &ReferenceImplementations::comm220ss);
  ReferenceImplementations.def("comm111ss",
                               &ReferenceImplementations::comm111ss);
  ReferenceImplementations.def("comm121ss",
                               &ReferenceImplementations::comm121ss);
  ReferenceImplementations.def("comm221ss",
                               &ReferenceImplementations::comm221ss);
  ReferenceImplementations.def("comm122ss",
                               &ReferenceImplementations::comm122ss);
  ReferenceImplementations.def("comm222_pp_hh_221ss",
                               &ReferenceImplementations::comm222_pp_hh_221ss);
  ReferenceImplementations.def("comm222_pp_hhss",
                               &ReferenceImplementations::comm222_pp_hhss);
  ReferenceImplementations.def("comm222_phss",
                               &ReferenceImplementations::comm222_phss);

  ReferenceImplementations.def("comm111st",
                               &ReferenceImplementations::comm111st);
  ReferenceImplementations.def("comm121st",
                               &ReferenceImplementations::comm121st);
  ReferenceImplementations.def("comm122st",
                               &ReferenceImplementations::comm122st);
  ReferenceImplementations.def("comm221st",
                               &ReferenceImplementations::comm221st);
  ReferenceImplementations.def("comm222_pp_hhst",
                               &ReferenceImplementations::comm222_pp_hhst);
  ReferenceImplementations.def("comm222_phst",
                               &ReferenceImplementations::comm222_phst);

  //
  ReferenceImplementations.def("comm331ss",
                               &ReferenceImplementations::comm331ss);
  ReferenceImplementations.def("comm223ss",
                               &ReferenceImplementations::comm223ss);
  ReferenceImplementations.def("comm231ss",
                               &ReferenceImplementations::comm231ss);
  ReferenceImplementations.def("comm232ss",
                               &ReferenceImplementations::comm232ss);
  ReferenceImplementations.def("comm133ss",
                               &ReferenceImplementations::comm133ss);
  ReferenceImplementations.def("comm132ss",
                               &ReferenceImplementations::comm132ss);
  ReferenceImplementations.def("comm332_ppph_hhhpss",
                               &ReferenceImplementations::comm332_ppph_hhhpss);
  ReferenceImplementations.def("comm332_pphhss",
                               &ReferenceImplementations::comm332_pphhss);
  ReferenceImplementations.def("comm233_pp_hhss",
                               &ReferenceImplementations::comm233_pp_hhss);
  ReferenceImplementations.def("comm233_phss",
                               &ReferenceImplementations::comm233_phss);
  ReferenceImplementations.def("comm333_ppp_hhhss",
                               &ReferenceImplementations::comm333_ppp_hhhss);
  ReferenceImplementations.def("comm333_pph_hhpss",
                               &ReferenceImplementations::comm333_pph_hhpss);

  //
  ReferenceImplementations.def("diagram_CIa",
                               &ReferenceImplementations::diagram_CIa);
  ReferenceImplementations.def("diagram_CIb",
                               &ReferenceImplementations::diagram_CIb);
  ReferenceImplementations.def("diagram_CIIa",
                               &ReferenceImplementations::diagram_CIIa);
  ReferenceImplementations.def("diagram_CIIb",
                               &ReferenceImplementations::diagram_CIIb);
  ReferenceImplementations.def("diagram_CIIc",
                               &ReferenceImplementations::diagram_CIIc);
  ReferenceImplementations.def("diagram_CIId",
                               &ReferenceImplementations::diagram_CIId);
  ReferenceImplementations.def("diagram_CIIIa",
                               &ReferenceImplementations::diagram_CIIIa);
  ReferenceImplementations.def("diagram_CIIIb",
                               &ReferenceImplementations::diagram_CIIIb);
  ReferenceImplementations.def("diagram_DIa",
                               &ReferenceImplementations::diagram_DIa);
  ReferenceImplementations.def("diagram_DIb",
                               &ReferenceImplementations::diagram_DIb);
  ReferenceImplementations.def("diagram_DIVa",
                               &ReferenceImplementations::diagram_DIVa);
  ReferenceImplementations.def("diagram_DIVb",
                               &ReferenceImplementations::diagram_DIVb);
  ReferenceImplementations.def(
      "diagram_DIVb_intermediate",
      &ReferenceImplementations::diagram_DIVb_intermediate);
  ReferenceImplementations.def(
      "comm223_231_BruteForce",
      &ReferenceImplementations::comm223_231_BruteForce);
  ReferenceImplementations.def(
      "comm223_232_BruteForce",
      &ReferenceImplementations::comm223_232_BruteForce);
  ReferenceImplementations.def("comm223_231",
                               &ReferenceImplementations::comm223_231);
  ReferenceImplementations.def("comm223_232",
                               &ReferenceImplementations::comm223_232);
  ReferenceImplementations.def("comm223_231_tts",
                               &ReferenceImplementations::comm223_231_tts);
  ReferenceImplementations.def("comm223_231_tts_fI",
                               &ReferenceImplementations::comm223_231_tts_fI);
  ReferenceImplementations.def("comm223_231_tts_fII",
                               &ReferenceImplementations::comm223_231_tts_fII);
  ReferenceImplementations.def("comm223_231_tts_fIIIa",
                               &ReferenceImplementations::comm223_231_tts_fIIIa);
  ReferenceImplementations.def("comm223_231_tts_fIIIb",
                               &ReferenceImplementations::comm223_231_tts_fIIIb);
  ReferenceImplementations.def("comm223_232_tts",
                               &ReferenceImplementations::comm223_232_tts);
  ReferenceImplementations.def("comm223_232_tts_GI",
                               &ReferenceImplementations::comm223_232_tts_GI);
  ReferenceImplementations.def("comm223_232_tts_GII",
                               &ReferenceImplementations::comm223_232_tts_GII);
  ReferenceImplementations.def("comm223_232_tts_GIIIb",
                               &ReferenceImplementations::comm223_232_tts_GIIIb);
  ReferenceImplementations.def("comm223_232_tts_GIIIc",
                               &ReferenceImplementations::comm223_232_tts_GIIIc,
                               py::arg("Eta"), py::arg("Gamma"), py::arg("Z"),
                               py::arg("which_term") = 0);
  ReferenceImplementations.def("comm223_232_tts_GIIIc_term1",
                               &ReferenceImplementations::comm223_232_tts_GIIIc_term1);
  ReferenceImplementations.def("comm223_232_tts_GIIIc_term2",
                               &ReferenceImplementations::comm223_232_tts_GIIIc_term2);
  ReferenceImplementations.def(
      "comm223_232_tts_GIIIc_tensor_red",
      &ReferenceImplementations::comm223_232_tts_GIIIc_tensor_red,
      py::arg("Eta"), py::arg("Gamma"), py::arg("Z"),
      py::arg("which_term") = 0);
  ReferenceImplementations.def("comm223_232_tts_GIVa",
                               &ReferenceImplementations::comm223_232_tts_GIVa);
  ReferenceImplementations.def("comm223_232_tts_GIVb",
                               &ReferenceImplementations::comm223_232_tts_GIVb,
                               py::arg("Eta"), py::arg("Gamma"), py::arg("Z"),
                               py::arg("which_term") = 0);
  ReferenceImplementations.def("comm223_232_tts_GIVc",
                               &ReferenceImplementations::comm223_232_tts_GIVc);
  ReferenceImplementations.def(
      "comm223_231_tts_BruteForce",
      &ReferenceImplementations::comm223_231_tts_BruteForce);
  ReferenceImplementations.def(
      "comm223_232_tts_BruteForce",
      &ReferenceImplementations::comm223_232_tts_BruteForce);
     ReferenceImplementations.def("comm223_132",
                                                                            &ReferenceImplementations::comm223_132);
     ReferenceImplementations.def("comm223_132_ladder",
                                                                            &ReferenceImplementations::comm223_132_ladder);
     ReferenceImplementations.def("comm223_132_cross",
                                                                            &ReferenceImplementations::comm223_132_cross);
     ReferenceImplementations.def("comm223_132_onebody",
                                                                            &ReferenceImplementations::comm223_132_onebody);
  ReferenceImplementations.def("evc_z1_mscheme",
                               &ReferenceImplementations::evc_z1_mscheme);
  ReferenceImplementations.def("evc_z2_mscheme",
                               &ReferenceImplementations::evc_z2_mscheme);
  ReferenceImplementations.def("evc_z0_mscheme",
                               &ReferenceImplementations::evc_z0_mscheme);

  ReferenceImplementations.def("comm331st",
                               &ReferenceImplementations::comm331st);
  ReferenceImplementations.def("comm223st",
                               &ReferenceImplementations::comm223st);
  ReferenceImplementations.def("comm231st",
                               &ReferenceImplementations::comm231st);
  ReferenceImplementations.def("comm232st",
                               &ReferenceImplementations::comm232st);
  ReferenceImplementations.def("comm232st_amc",
                               &ReferenceImplementations::comm232st_amc);
  ReferenceImplementations.def("comm232st_amc_eq1",
                               &ReferenceImplementations::comm232st_amc_eq1);
  ReferenceImplementations.def("comm133st",
                               &ReferenceImplementations::comm133st);
  ReferenceImplementations.def("comm132st",
                               &ReferenceImplementations::comm132st);
  ReferenceImplementations.def("comm231tts",
                               &ReferenceImplementations::comm231tts);
  ReferenceImplementations.def("comm110tts",
                               &ReferenceImplementations::comm110tts);
  ReferenceImplementations.def("comm220tts",
                               &ReferenceImplementations::comm220tts);
  ReferenceImplementations.def("comm111tts",
                               &ReferenceImplementations::comm111tts);
  ReferenceImplementations.def("comm121tts",
                               &ReferenceImplementations::comm121tts);
  ReferenceImplementations.def("comm122tts",
                               &ReferenceImplementations::comm122tts);
  ReferenceImplementations.def("comm221tts",
                               &ReferenceImplementations::comm221tts);
  ReferenceImplementations.def("comm222_pp_hhtts",
                               &ReferenceImplementations::comm222_pp_hhtts);
  ReferenceImplementations.def("comm222_phtts",
                               &ReferenceImplementations::comm222_phtts);
  ReferenceImplementations.def("comm132tts",
                               &ReferenceImplementations::comm132tts);
  ReferenceImplementations.def("comm232tts",
                               &ReferenceImplementations::comm232tts);
  ReferenceImplementations.def("comm232tts_bare",
                               &ReferenceImplementations::comm232tts_bare);
  ReferenceImplementations.def("comm223tts",
                               &ReferenceImplementations::comm223tts);
  ReferenceImplementations.def("comm332_ppph_hhhpst",
                               &ReferenceImplementations::comm332_ppph_hhhpst);
  ReferenceImplementations.def("comm332_pphhst",
                               &ReferenceImplementations::comm332_pphhst);
  ReferenceImplementations.def("comm233_pp_hhst",
                               &ReferenceImplementations::comm233_pp_hhst);
  ReferenceImplementations.def("comm233_phst",
                               &ReferenceImplementations::comm233_phst);
  ReferenceImplementations.def("comm333_ppp_hhhst",
                               &ReferenceImplementations::comm333_ppp_hhhst);
  ReferenceImplementations.def("comm333_pph_hhpst",
                               &ReferenceImplementations::comm333_pph_hhpst);
     ReferenceImplementations.def("TriplesGuess",
                                                                            &ReferenceImplementations::TriplesGuess);

  // EOM::ArnoldiResult — returned by ArnoldiSolve and Run()
  py::class_<EOM::ArnoldiResult>(m, "ArnoldiResult")
      .def_property_readonly("energies",
           [](const EOM::ArnoldiResult &r) {
             std::vector<double> v(r.energies.begin(), r.energies.end());
             return v;
                          })
               .def_property_readonly("eigvecs",
                          [](const EOM::ArnoldiResult &r) {
                               std::vector<std::vector<double>> out(
                                         r.eigvecs.n_rows,
                                         std::vector<double>(r.eigvecs.n_cols, 0.0));
                               for (size_t i = 0; i < r.eigvecs.n_rows; ++i)
                                    for (size_t j = 0; j < r.eigvecs.n_cols; ++j)
                                         out[i][j] = r.eigvecs(i, j);
                               return out;
                          })
               .def_property_readonly("hall",
                          [](const EOM::ArnoldiResult &r) {
                               std::vector<std::vector<double>> out(
                                         r.hall.n_rows,
                                         std::vector<double>(r.hall.n_cols, 0.0));
                               for (size_t i = 0; i < r.hall.n_rows; ++i)
                                    for (size_t j = 0; j < r.hall.n_cols; ++j)
                                         out[i][j] = r.hall(i, j);
                               return out;
                          })
               .def_property_readonly("ritz",
                          [](const EOM::ArnoldiResult &r) {
                               return r.ritz;
                          })
               .def_property_readonly("residuals",
                          [](const EOM::ArnoldiResult &r) {
                               std::vector<double> v(r.residuals.n_elem);
                               for (size_t i = 0; i < r.residuals.n_elem; ++i)
                                    v[i] = r.residuals(i);
                               return v;
                          })
               .def_readonly("max_ortho", &EOM::ArnoldiResult::max_ortho)
               .def_readonly("steps", &EOM::ArnoldiResult::steps)
               .def_readonly("converged", &EOM::ArnoldiResult::converged)
               .def_readonly("stop_reason", &EOM::ArnoldiResult::stop_reason);

     py::class_<EOM::ArnoldiTraceDiffResult>(m, "ArnoldiTraceDiffResult")
               .def_readonly("found", &EOM::ArnoldiTraceDiffResult::found)
               .def_readonly("step", &EOM::ArnoldiTraceDiffResult::step)
               .def_readonly("i", &EOM::ArnoldiTraceDiffResult::i)
               .def_readonly("j", &EOM::ArnoldiTraceDiffResult::j)
               .def_readonly("hall_new", &EOM::ArnoldiTraceDiffResult::hall_new)
               .def_readonly("hall_old", &EOM::ArnoldiTraceDiffResult::hall_old)
               .def_readonly("delta_hall", &EOM::ArnoldiTraceDiffResult::delta_hall)
               .def_readonly("h1_sym_new", &EOM::ArnoldiTraceDiffResult::h1_sym_new)
               .def_readonly("h1_sym_old", &EOM::ArnoldiTraceDiffResult::h1_sym_old)
               .def_readonly("h2_cross_new", &EOM::ArnoldiTraceDiffResult::h2_cross_new)
               .def_readonly("h2_cross_old", &EOM::ArnoldiTraceDiffResult::h2_cross_old)
               .def_readonly("h2_sym_new", &EOM::ArnoldiTraceDiffResult::h2_sym_new)
               .def_readonly("h2_sym_old", &EOM::ArnoldiTraceDiffResult::h2_sym_old)
               .def_readonly("h2_diag_i_new", &EOM::ArnoldiTraceDiffResult::h2_diag_i_new)
               .def_readonly("h2_diag_j_new", &EOM::ArnoldiTraceDiffResult::h2_diag_j_new)
               .def_readonly("h2_diag_i_old", &EOM::ArnoldiTraceDiffResult::h2_diag_i_old)
               .def_readonly("h2_diag_j_old", &EOM::ArnoldiTraceDiffResult::h2_diag_j_old)
               .def_readonly("max_abs_diff", &EOM::ArnoldiTraceDiffResult::max_abs_diff)
               .def_readonly("max_rel_diff", &EOM::ArnoldiTraceDiffResult::max_rel_diff);

  // EOM::RunResult — returned by EOM.Run()
  py::class_<EOM::RunResult>(m, "RunResult")
      .def_readonly("eref",    &EOM::RunResult::eref)
      .def_readonly("arnoldi", &EOM::RunResult::arnoldi)
      .def_readonly("Q_orbit", &EOM::RunResult::Q_orbit,
           "PA/PR: attached particle or removed hole orbit index (-1 if N/A)")
      .def_readonly("spe", &EOM::RunResult::spe,
           "Hs OneBody(Q,Q) in the PA/PR channel")
      .def_readonly("rayleigh_1h", &EOM::RunResult::rayleigh_1h,
           "PR: 1h Rayleigh quotient on -[H,a†] (≈ -SPE)");

  py::enum_<EOM::SREOMMode>(m, "SREOMMode")
      .value("Excitation", EOM::SREOMMode::Excitation)
      .value("ParticleAttached", EOM::SREOMMode::ParticleAttached)
      .value("ParticleRemoved", EOM::SREOMMode::ParticleRemoved)
      .export_values();

  py::class_<EOM>(m, "EOM")
      // constructors
      .def(py::init<Operator &, Operator &, int, int, int>(),
           py::arg("Hs"), py::arg("rdm"), py::arg("J2"), py::arg("parity"), py::arg("itz"))
      .def(py::init<Operator &, const std::string &, int, int, int>(),
           py::arg("Hs"), py::arg("tdm_file"), py::arg("J2"), py::arg("parity"), py::arg("itz"))
      .def(py::init<Operator &, int, int, int>(),
           py::arg("Hs"), py::arg("J2"), py::arg("parity"), py::arg("itz"))
      .def(py::init<Operator &, int, int, int, EOM::SREOMMode>(),
           py::arg("Hs"), py::arg("J2"), py::arg("parity"), py::arg("itz"),
           py::arg("sr_mode"),
           "Single-reference EOM; sr_mode = Excitation | ParticleAttached | ParticleRemoved")
      // high-level entry point: init + solve in one call
      .def("Run", &EOM::Run,
           py::arg("max_iter") = 200, py::arg("state_want") = 6,
           "Solve: dispatches on is_multiref and sr_mode set at construction or via SetSREOMMode")
      .def("RunSR", &EOM::RunSR, py::arg("max_iter") = 200, py::arg("state_want") = 6,
           "SR same-A excitation EOM (1p1h ⊕ 2p2h)")
      .def("RunPA", &EOM::RunPA, py::arg("max_iter") = 200, py::arg("state_want") = 6,
           "SR 1-particle-attached EOM (A+1); ω = E(A+1) - E(A)")
      .def("RunPR", &EOM::RunPR, py::arg("max_iter") = 200, py::arg("state_want") = 6,
           "SR 1-particle-removed EOM (A-1); ω = E(A-1) - E(A)")
      .def("SetSREOMMode", &EOM::SetSREOMMode, py::arg("mode"),
           "SR only: Excitation (default), ParticleAttached, ParticleRemoved")
      .def("GetSREOMMode", &EOM::GetSREOMMode)
      // MR setup — must be called before operator-level methods in MR mode
      .def("ComputeNorm", &EOM::ComputeNorm, py::arg("Op1"), py::arg("Op2"))
      .def("FlattenOperator", [](const EOM &self, Operator &Op) {
            arma::vec v = self.FlattenOperator(Op);
            return py::array_t<double>({(py::ssize_t)v.n_elem}, v.memptr());
          }, py::arg("Op"))
      .def("UnflattenOperator", [](const EOM &self, Operator &Op,
                                   py::array_t<double, py::array::c_style | py::array::forcecast> arr) {
            py::buffer_info info = arr.request();
            if (info.ndim != 1)
              throw std::runtime_error("UnflattenOperator: need 1-D array");
            arma::vec v(info.shape[0], arma::fill::zeros);
            const double *p = static_cast<double *>(info.ptr);
            for (py::ssize_t i = 0; i < info.shape[0]; ++i)
              v((arma::uword)i) = p[i];
            self.UnflattenOperator(Op, v);
          }, py::arg("Op"), py::arg("v"))
      .def("ConstructConfigs",       &EOM::ConstructConfigs)
      .def("ConstructNormMatrix",    &EOM::ConstructNormMatrix)
      .def("ConstructProjectMatrix", &EOM::ConstructProjectMatrix)
      .def("SetArnoldiUseProjection", &EOM::SetArnoldiUseProjection,
           py::arg("use_projection"))
      .def("SetArnoldiCheckExpectation", &EOM::SetArnoldiCheckExpectation,
           py::arg("check_expectation"))
      .def("SetArnoldiPrintTiming", &EOM::SetArnoldiPrintTiming,
           py::arg("print_timing"))
      .def("SetArnoldiMonitorOrtho", &EOM::SetArnoldiMonitorOrtho,
           py::arg("monitor"))
      .def("SetArnoldiUseH3", &EOM::SetArnoldiUseH3, py::arg("use_h3"))
      .def("SetUseRdm3", &EOM::SetUseRdm3, py::arg("use"),
           "Include valence ρ₃ (TRBTD) in N / overlap; default false")
      .def("GetUseRdm3", &EOM::GetUseRdm3)
      .def("SetIncludeConfigs", &EOM::SetIncludeConfigs,
           py::arg("qv"), py::arg("ph"), py::arg("ppvv"), py::arg("pphv"),
           py::arg("pphh"),
           "Select EOM config blocks (qv,ph,ppvv,pphv,pphh); call before ConstructConfigs")
      .def("PrintIncludeConfigs", &EOM::PrintIncludeConfigs)
      .def("SetArnoldiPreferPositive", &EOM::SetArnoldiPreferPositive,
           py::arg("prefer"),
           "Prefer ΔE >= soft_floor when packing returned roots")
      .def("SetArnoldiSoftFloor", &EOM::SetArnoldiSoftFloor, py::arg("floor_mev"),
           "Soft-mode floor (MeV) for prefer-positive selection")
      .def("SetArnoldiEnergyTol", &EOM::SetArnoldiEnergyTol, py::arg("tol"))
      .def("SetArnoldiResidTol", &EOM::SetArnoldiResidTol, py::arg("tol"))
      .def("SetArnoldiOrthoWarn", &EOM::SetArnoldiOrthoWarn, py::arg("tol"))
      .def("SetArnoldiOrthoFail", &EOM::SetArnoldiOrthoFail, py::arg("tol"))
      .def("SetArnoldiStableWindow", &EOM::SetArnoldiStableWindow, py::arg("n"))
      .def("SetReferenceEnergyShift", &EOM::SetReferenceEnergyShift,
           py::arg("energy"))
      .def("ClearReferenceEnergyShift", &EOM::ClearReferenceEnergyShift)
      .def("GetReferenceEnergyShift", &EOM::GetReferenceEnergyShift)
      .def("PrepareHamiltonianForArnoldi", &EOM::PrepareHamiltonianForArnoldi)
      // operator-level building blocks used in Python EOM loops
      .def("force_decouple",            &EOM::force_decouple,
           py::arg("H"))
      .def("EraseValence",             &EOM::EraseValence,
           py::arg("H"),
           "DISABLED no-op: dropping vv + E_val*N is incorrect for MR-EOM")
      .def("EraseQspace",              &EOM::EraseQspace,
           py::arg("H"))
      .def("ProjectOprator",            &EOM::ProjectOprator,
           py::arg("Qin"))
      .def("GetVSEOM_Overlap_single",   &EOM::GetVSEOM_Overlap_single,
           py::arg("H1"), py::arg("H2"))
      .def("GetVSEOM_Overlap_multiref", &EOM::GetVSEOM_Overlap_multiref,
           py::arg("H"))
      .def("NormMultiref",             &EOM::NormMultiref,
           py::arg("T1"), py::arg("T2"))
      .def("Norm3Multiref",            &EOM::Norm3Multiref,
           py::arg("t1"), py::arg("t2"), py::arg("haml"))
      .def("HtcMultiref",              &EOM::HtcMultiref,
           py::arg("haml"), py::arg("chi"))
      .def("DcomMultiref",             &EOM::DcomMultiref,
           py::arg("haml"), py::arg("chi"))
      .def("ExpectationValue",         &EOM::ExpectationValue,
           py::arg("Psi"),
           "Delta E from full decoupled Hs (vv kept); no E_val shift")
      .def("ExpectationValueFull",     &EOM::ExpectationValueFull,
           py::arg("Psi"), py::arg("use_rdm_norm") = false,
           "Delta E from Hs_full (vv kept); no E_val shift. "
           "Dcom already includes the S vs chi 1/2.")
      .def("HaveHsFull",               &EOM::HaveHsFull)
      .def("GetHsFull",                &EOM::GetHsFull,
           py::return_value_policy::reference_internal)
      .def("ArnoldiSolve",     &EOM::ArnoldiSolve,
           py::arg("vi"), py::arg("max_iter"), py::arg("state_want"))
      .def("BuildCanonicalTransform", &EOM::BuildCanonicalTransform,
           py::arg("eps") = 1e-8,
           "Canonical ortho X=U s^{-1/2} of Nkernel; returns M retained")
      .def("SolveGEPCanonical", &EOM::SolveGEPCanonical,
           py::arg("state_want"), py::arg("eps") = 1e-8,
           py::arg("n_h2_subspace") = 50,
           "Dense GEP: canon. ortho; if UseH3, H2 then H2+H3 in n_h2_subspace")
      .def("ArnoldiSolveH2",   &EOM::ArnoldiSolveH2,
           py::arg("vi"), py::arg("max_iter"), py::arg("state_want"),
           "Arnoldi on full decoupled H (no H_PP/H_QQ split, no E_val shift); "
           "Hall = (Htc_ab+Htc_ba)/2 + Dcom (default; matches ExpectationValue)")
      .def("ArnoldiSolve_old", &EOM::ArnoldiSolve_old,
           py::arg("vi"), py::arg("max_iter"), py::arg("state_want"))
      .def("CompareArnoldiHallBuild", &EOM::CompareArnoldiHallBuild,
           py::arg("vi"), py::arg("max_iter"), py::arg("tol") = 1e-10)
      .def("GetVSEOM_ladder_single",    &EOM::GetVSEOM_ladder_single,
           py::arg("H"), py::arg("herm"))
      .def("GetVSEOM_ladder_multiref",  &EOM::GetVSEOM_ladder_multiref,
           py::arg("H"), py::arg("herm"))
      .def("ReadTdm", &EOM::ReadTdm,
           py::arg("tdm_file"))
      .def("WriteTdm", &EOM::WriteTdm,
           py::arg("op"), py::arg("filename"))
      .def_readonly("ppvv_start", &EOM::ppvv_start)
      .def_readonly("ppvv_end",   &EOM::ppvv_end)
      .def_readonly("ppvv_dim",   &EOM::ppvv_dim)
      .def_readonly("pphv_start", &EOM::pphv_start)
      .def_readonly("pphv_end",   &EOM::pphv_end)
      .def_readonly("pphv_dim",   &EOM::pphv_dim)
      .def_readonly("pphh_start", &EOM::pphh_start)
      .def_readonly("pphh_end",   &EOM::pphh_end)
      .def_readonly("eom_dims",   &EOM::eom_dims)
      .def_property_readonly("eom_confs", [](const EOM &self) {
            std::vector<std::array<size_t,4>> out;
            for (auto &c : self.eom_confs)
              out.push_back({(size_t)c[0], (size_t)c[1], (size_t)c[2], (size_t)c[3]});
            return out;
          })
      .def("ThreeBody_Diagram", &EOM::ThreeBody_Diagram,
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"),
           py::arg("f"), py::arg("g"), py::arg("j0"), py::arg("j2"))
      .def("ThreeBody_Diagram_Entries", &EOM::ThreeBody_Diagram_Entries,
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), py::arg("e"),
           py::arg("f"), py::arg("g"), py::arg("j0"), py::arg("j2"))
      .def("RdmThreeBody_J", &EOM::RdmThreeBody_J,
           py::arg("Jab"), py::arg("a"), py::arg("b"), py::arg("c"),
           py::arg("Jde"), py::arg("d"), py::arg("e"), py::arg("f"),
           py::arg("twoJ"))
      .def("ShowModel", &EOM::ShowModel)
      .def_readonly("rdm", &EOM::rdm)
      .def("PrintConfigs", &EOM::PrintConfigs)
      .def_property_readonly("Nkernel", [](const EOM &self) {
          // Convert arma::sp_mat to dense numpy array
          arma::mat dense = arma::mat(self.Nkernel);
          return py::array_t<double>(
              {(py::ssize_t)dense.n_rows, (py::ssize_t)dense.n_cols},
              dense.memptr());
      });

  py::class_<RPA>(m, "RPA")
      .def(py::init<Operator &>())
      .def("ConstructAMatrix", &RPA::ConstructAMatrix, py::arg("J"),
           py::arg("parity"), py::arg("Tz"), py::arg("Isovector"))
      .def("ConstructBMatrix", &RPA::ConstructBMatrix, py::arg("J"),
           py::arg("parity"), py::arg("Tz"), py::arg("Isovector"))
      .def("SolveCP", &RPA::SolveCP)
      .def("SolveTDA", &RPA::SolveTDA)
      .def("SolveRPA", &RPA::SolveRPA)
      .def("TransitionToGroundState", &RPA::TransitionToGroundState,
           py::arg("OpIn"), py::arg("mu"))
      .def("PVCouplingEffectiveCharge", &RPA::PVCouplingEffectiveCharge,
           py::arg("OpIn"), py::arg("k"), py::arg("l"))
      .def("GetEnergies",
           [](RPA &self) {
             arma::vec vals = self.GetEnergies();
             std::vector<double> vvec;
             for (auto &v : vals) {
               vvec.push_back(v);
             };
             return vvec;
           })
      .def("GetX",
           [](RPA &self, size_t i) {
             arma::vec vals = self.GetX(i);
             std::vector<double> vvec;
             for (auto &v : vals) {
               vvec.push_back(v);
             };
             return vvec;
           })
      .def("GetY",
           [](RPA &self, size_t i) {
             arma::vec vals = self.GetY(i);
             std::vector<double> vvec;
             for (auto &v : vals) {
               vvec.push_back(v);
             };
             return vvec;
           })
      .def("PrintA", [](RPA &self) { std::cout << self.A << std::endl; })
      .def("PrintB", [](RPA &self) { std::cout << self.B << std::endl; })
      .def("GetEgs", &RPA::GetEgs);

  py::class_<EVC> evc(m, "EVC");
  py::class_<EVC::ClusterAmplitudes>(evc, "ClusterAmplitudes")
      .def_readwrite("z0", &EVC::ClusterAmplitudes::z0)
      .def_readwrite("z", &EVC::ClusterAmplitudes::z);
  py::class_<EVC::RHS>(evc, "RHS")
      .def_readwrite("dz0", &EVC::RHS::dz0)
      .def_readwrite("dz", &EVC::RHS::dz);
  evc.def(py::init<ModelSpace &>())
      .def("SetEulerSteps", &EVC::SetEulerSteps)
      .def("SetUseRK4", &EVC::SetUseRK4)
      .def("ExtractExcitationPart", &EVC::ExtractExcitationPart)
      .def("ExtractTFromOmega", &EVC::ExtractTFromOmega)
      .def("BuildZ0RHS", &EVC::BuildZ0RHS)
      .def("BuildZ1RHS", &EVC::BuildZ1RHS)
      .def("BuildZ2RHS", &EVC::BuildZ2RHS)
      .def("BuildRHS", &EVC::BuildRHS)
      .def("Solve", &EVC::Solve, py::arg("t"), py::arg("lambda_max") = 1.0)
      .def("SolveFromOmega", &EVC::SolveFromOmega, py::arg("omega"), py::arg("lambda_max") = 1.0)
      .def("SolveEuler", &EVC::SolveEuler, py::arg("t"), py::arg("tdagger"), py::arg("lambda_max") = 1.0)
      .def("CoupledClusterEnergy", &EVC::CoupledClusterEnergy)
      .def("HamiltonianKernel", &EVC::HamiltonianKernel)
      .def("NormKernel", &EVC::NormKernel)
      .def("EvaluateKernels", [](EVC &self, const Operator &h0, const Operator &omega_i, const Operator &omega_k)
           {
             const auto kn = self.EvaluateKernels(h0, omega_i, omega_k);
             return py::make_tuple(kn[0], kn[1]);
           });

  py::class_<UnitTest>(m, "UnitTest")
      //      .def(py::init<>())
      .def(py::init<ModelSpace &>())
      .def("SetRandomSeed", &UnitTest::SetRandomSeed)
      .def("RandomDaggerOp", &UnitTest::RandomDaggerOp, py::arg("modelspace"),
           py::arg("Q"))
      .def("RandomOp", &UnitTest::RandomOp, py::arg("modelspace"),
           py::arg("jrank"), py::arg("tz"), py::arg("parity"),
           py::arg("particle_rank"), py::arg("hermitian"))
      .def("TestCommutators", &UnitTest::TestCommutators)
      .def("TestCommutators_Tensor", &UnitTest::TestCommutators_Tensor)
      .def("TestCommutators_IsospinChanging",
           &UnitTest::TestCommutators_IsospinChanging)
      .def("TestCommutators_ParityChanging",
           &UnitTest::TestCommutators_ParityChanging)
      .def("TestCommutators3", &UnitTest::TestCommutators3)
      .def("TestNormalOrdering", &UnitTest::TestNormalOrdering)
      .def("TestDaggerCommutators", &UnitTest::TestDaggerCommutators)
      .def("TestDaggerCommutatorsAlln", &UnitTest::TestDaggerCommutatorsAlln)
      .def("Test_comm211sd", &UnitTest::Test_comm211sd)
      .def("Test_comm231sd", &UnitTest::Test_comm231sd)
      .def("Test_comm431sd", &UnitTest::Test_comm431sd)
      .def("Test_comm413sd", &UnitTest::Test_comm413sd)
      .def("Test_comm233sd", &UnitTest::Test_comm233sd)
      .def("Test_comm433_pp_hh_sd", &UnitTest::Test_comm433_pp_hh_sd)
      .def("Test_comm433sd_ph", &UnitTest::Test_comm433sd_ph)
      .def("Test3BodyAntisymmetry", &UnitTest::Test3BodyAntisymmetry)
      .def("Test3BodyHermiticity", &UnitTest::Test3BodyHermiticity)
      .def("TestRPAEffectiveCharge", &UnitTest::TestRPAEffectiveCharge,
           py::arg("H"), py::arg("OpIn"), py::arg("k"), py::arg("l"))
      .def("SanityCheck", &UnitTest::SanityCheck)
      .def("TestFactorizedDoubleCommutators",
           &UnitTest::TestFactorizedDoubleCommutators)
      .def("TestPerturbativeTriples", &UnitTest::TestPerturbativeTriples)
      .def("Test_evc_rhs_ccsd", &UnitTest::Test_evc_rhs_ccsd)
      .def("Test_evc_z1_jscheme", &UnitTest::Test_evc_z1_jscheme)
      .def("Test_evc_z2_jscheme", &UnitTest::Test_evc_z2_jscheme)
      .def("Test_evc_z0_jscheme", &UnitTest::Test_evc_z0_jscheme)
      .def("Test_evc_ode", &UnitTest::Test_evc_ode)
      .def("Test_evc_kernels", &UnitTest::Test_evc_kernels)
      .def("Test_comm110ss", &UnitTest::Test_comm110ss)
      .def("Test_comm220ss", &UnitTest::Test_comm220ss)
      .def("Test_comm111ss", &UnitTest::Test_comm111ss)
      .def("Test_comm121ss", &UnitTest::Test_comm121ss)
      .def("Test_comm221ss", &UnitTest::Test_comm221ss)
      .def("Test_comm122ss", &UnitTest::Test_comm122ss)
      .def("Test_comm222_pp_hhss", &UnitTest::Test_comm222_pp_hhss)
      .def("Test_comm222_phss", &UnitTest::Test_comm222_phss)
      .def("Test_comm222_pp_hh_221ss", &UnitTest::Test_comm222_pp_hh_221ss)

      .def("Test_comm111st", &UnitTest::Test_comm111st)
      .def("Test_comm121st", &UnitTest::Test_comm121st)
      .def("Test_comm221st", &UnitTest::Test_comm221st)
      .def("Test_comm122st", &UnitTest::Test_comm122st)
      .def("Test_comm222_pp_hhst", &UnitTest::Test_comm222_pp_hhst)
      .def("Test_comm222_phst", &UnitTest::Test_comm222_phst)

      .def("Test_comm330ss", &UnitTest::Test_comm330ss)
      .def("Test_comm331ss", &UnitTest::Test_comm331ss)
      .def("Test_comm231ss", &UnitTest::Test_comm231ss)
      .def("Test_comm132ss", &UnitTest::Test_comm132ss)
      .def("Test_comm232ss", &UnitTest::Test_comm232ss)
      .def("Test_comm223ss", &UnitTest::Test_comm223ss)
      .def("Test_comm133ss", &UnitTest::Test_comm133ss)
      .def("Test_comm332_ppph_hhhpss", &UnitTest::Test_comm332_ppph_hhhpss)
      .def("Test_comm332_pphhss", &UnitTest::Test_comm332_pphhss)
      .def("Test_comm233_pp_hhss", &UnitTest::Test_comm233_pp_hhss)
      .def("Test_comm233_phss", &UnitTest::Test_comm233_phss)
      .def("Test_comm333_ppp_hhhss", &UnitTest::Test_comm333_ppp_hhhss)
      .def("Test_comm333_pph_hhpss", &UnitTest::Test_comm333_pph_hhpss)

      .def("Mscheme_Test_comm110ss", &UnitTest::Mscheme_Test_comm110ss)
      .def("Mscheme_Test_comm220ss", &UnitTest::Mscheme_Test_comm220ss)
      .def("Mscheme_Test_comm111ss", &UnitTest::Mscheme_Test_comm111ss)
      .def("Mscheme_Test_comm121ss", &UnitTest::Mscheme_Test_comm121ss)
      .def("Mscheme_Test_comm221ss", &UnitTest::Mscheme_Test_comm221ss)
      .def("Mscheme_Test_comm122ss", &UnitTest::Mscheme_Test_comm122ss)
      .def("Mscheme_Test_comm222_pp_hhss",
           &UnitTest::Mscheme_Test_comm222_pp_hhss)
      .def("Mscheme_Test_comm222_phss", &UnitTest::Mscheme_Test_comm222_phss)
      //
      //      .def("Mscheme_Test_comm222_pp_hh_221ss",
      //      &UnitTest::Mscheme_Test_comm222_pp_hh_221ss)
      ///
      .def("Mscheme_Test_comm330ss", &UnitTest::Mscheme_Test_comm330ss)
      .def("Mscheme_Test_comm331ss", &UnitTest::Mscheme_Test_comm331ss)
      .def("Mscheme_Test_comm231ss", &UnitTest::Mscheme_Test_comm231ss)
      .def("Mscheme_Test_comm132ss", &UnitTest::Mscheme_Test_comm132ss)
      .def("Mscheme_Test_comm232ss", &UnitTest::Mscheme_Test_comm232ss)
      .def("Mscheme_Test_comm223ss", &UnitTest::Mscheme_Test_comm223ss)
      .def("Mscheme_Test_comm133ss", &UnitTest::Mscheme_Test_comm133ss)
      .def("Mscheme_Test_comm332_ppph_hhhpss",
           &UnitTest::Mscheme_Test_comm332_ppph_hhhpss)
      .def("Mscheme_Test_comm332_pphhss",
           &UnitTest::Mscheme_Test_comm332_pphhss)
      .def("Mscheme_Test_comm233_pp_hhss",
           &UnitTest::Mscheme_Test_comm233_pp_hhss)
      .def("Mscheme_Test_comm233_phss", &UnitTest::Mscheme_Test_comm233_phss)
      .def("Mscheme_Test_comm333_ppp_hhhss",
           &UnitTest::Mscheme_Test_comm333_ppp_hhhss)
      .def("Mscheme_Test_comm333_pph_hhpss",
           &UnitTest::Mscheme_Test_comm333_pph_hhpss)
      //      .def("Test3BodySetGet",&UnitTest::Test3BodySetGet)

      // Tensor commutator with 3b
      .def("Mscheme_Test_comm331st", &UnitTest::Mscheme_Test_comm331st)
      .def("Mscheme_Test_comm223st", &UnitTest::Mscheme_Test_comm223st)
      .def("Mscheme_Test_comm231st", &UnitTest::Mscheme_Test_comm231st)
      .def("Mscheme_Test_comm232st", &UnitTest::Mscheme_Test_comm232st)
      .def("Mscheme_Test_comm232st_amc", &UnitTest::Mscheme_Test_comm232st_amc)
      .def("Test_comm232st_amc", &UnitTest::Test_comm232st_amc)
      .def("Mscheme_Test_comm133st", &UnitTest::Mscheme_Test_comm133st)
      .def("Mscheme_Test_comm132st", &UnitTest::Mscheme_Test_comm132st)
      .def("Mscheme_Test_comm231tts", &UnitTest::Mscheme_Test_comm231tts)
      .def("Mscheme_Test_comm110tts", &UnitTest::Mscheme_Test_comm110tts)
      .def("Mscheme_Test_comm220tts", &UnitTest::Mscheme_Test_comm220tts)
      .def("Mscheme_Test_comm111tts", &UnitTest::Mscheme_Test_comm111tts)
      .def("Mscheme_Test_comm121tts", &UnitTest::Mscheme_Test_comm121tts)
      .def("Mscheme_Test_comm122tts", &UnitTest::Mscheme_Test_comm122tts)
      .def("Mscheme_Test_comm221tts", &UnitTest::Mscheme_Test_comm221tts)
      .def("Mscheme_Test_comm222_pp_hhtts", &UnitTest::Mscheme_Test_comm222_pp_hhtts)
      .def("Mscheme_Test_comm222_phtts", &UnitTest::Mscheme_Test_comm222_phtts)
      .def("Mscheme_Test_comm132tts", &UnitTest::Mscheme_Test_comm132tts)
      .def("Mscheme_Test_comm232tts", &UnitTest::Mscheme_Test_comm232tts)
      .def("Mscheme_Test_comm223tts", &UnitTest::Mscheme_Test_comm223tts)
      .def("RME_Test_comm232tts_bare", &UnitTest::RME_Test_comm232tts_bare)
      .def("Test_comm231tts", &UnitTest::Test_comm231tts)
      .def("Test_comm110tts", &UnitTest::Test_comm110tts)
      .def("Test_comm220tts", &UnitTest::Test_comm220tts)
      .def("Test_comm111tts", &UnitTest::Test_comm111tts)
      .def("Test_comm121tts", &UnitTest::Test_comm121tts)
      .def("Test_comm122tts", &UnitTest::Test_comm122tts)
      .def("Test_comm221tts", &UnitTest::Test_comm221tts)
      .def("Test_comm222_pp_hhtts", &UnitTest::Test_comm222_pp_hhtts)
      .def("Test_comm222_phtts", &UnitTest::Test_comm222_phtts)
      .def("Test_comm132tts", &UnitTest::Test_comm132tts)
      .def("Test_comm232tts", &UnitTest::Test_comm232tts)
      .def("Test_comm223tts", &UnitTest::Test_comm223tts)
      .def("Test_scalar_tts_matches_ss", &UnitTest::Test_scalar_tts_matches_ss)
      .def("Mscheme_Test_comm332_ppph_hhhpst",
           &UnitTest::Mscheme_Test_comm332_ppph_hhhpst)
      .def("Mscheme_Test_comm332_pphhst",
           &UnitTest::Mscheme_Test_comm332_pphhst)
      .def("Mscheme_Test_comm233_pp_hhst",
           &UnitTest::Mscheme_Test_comm233_pp_hhst)
      .def("Mscheme_Test_comm233_phst", &UnitTest::Mscheme_Test_comm233_phst)
      .def("Mscheme_Test_comm233_phst", &UnitTest::Mscheme_Test_comm233_phst)
      .def("Mscheme_Test_comm333_ppp_hhhst",
           &UnitTest::Mscheme_Test_comm333_ppp_hhhst)
      .def("Mscheme_Test_comm333_pph_hhpst",
           &UnitTest::Mscheme_Test_comm333_pph_hhpst)

      .def("GetMschemeMatrixElement_1b", &UnitTest::GetMschemeMatrixElement_1b,
           py::arg("Op"), py::arg("a"), py::arg("ma"), py::arg("b"),
           py::arg("mb")) // Op, a,ma, b,mb...
      .def("GetMschemeMatrixElement_2b",
           &UnitTest::GetMschemeMatrixElement_2b) // Op, a,ma, b,mb...
      .def("GetMschemeMatrixElement_3b",
           &UnitTest::GetMschemeMatrixElement_3b) // Op, a,ma, b,mb...
      .def("Mscheme_comm223_232_GIVb", &UnitTest::Mscheme_comm223_232_GIVb,
           py::arg("Eta"), py::arg("Gamma"),
           py::arg("i"), py::arg("mi"), py::arg("j"), py::arg("mj"),
           py::arg("k"), py::arg("mk"), py::arg("l"), py::arg("ml"),
           py::arg("which_term") = 0)
      .def("Mscheme_chi_alpha", &UnitTest::Mscheme_chi_alpha)
      .def("Mscheme_chi_beta", &UnitTest::Mscheme_chi_beta)
      .def("Mscheme_chi_gamma", &UnitTest::Mscheme_chi_gamma)
      .def("Mscheme_chi_delta", &UnitTest::Mscheme_chi_delta)
      .def("Mscheme_chi_epsilon", &UnitTest::Mscheme_chi_epsilon)
      .def("Mscheme_chi_zeta", &UnitTest::Mscheme_chi_zeta)
      .def("Mscheme_chi_OmegaGamma", &UnitTest::Mscheme_chi_OmegaGamma)
      .def("Mscheme_chi_eta", &UnitTest::Mscheme_chi_eta)
      .def("Mscheme_chi_theta", &UnitTest::Mscheme_chi_theta)
      .def("Mscheme_chi_iota", &UnitTest::Mscheme_chi_iota)
      .def("Mscheme_chi_kappa", &UnitTest::Mscheme_chi_kappa)
      .def("Mscheme_chi_lambda", &UnitTest::Mscheme_chi_lambda)
      .def("Mscheme_fact_fI", &UnitTest::Mscheme_fact_fI)
      .def("Mscheme_fact_fII", &UnitTest::Mscheme_fact_fII)
      .def("Mscheme_fact_fIIIa", &UnitTest::Mscheme_fact_fIIIa)
      .def("Mscheme_fact_fIIIb", &UnitTest::Mscheme_fact_fIIIb)
      .def("Mscheme_fact_GI", &UnitTest::Mscheme_fact_GI)
      .def("Mscheme_fact_GII", &UnitTest::Mscheme_fact_GII)
      .def("Mscheme_fact_GIIIa", &UnitTest::Mscheme_fact_GIIIa)
      .def("Mscheme_fact_GIIIb", &UnitTest::Mscheme_fact_GIIIb)
      .def("Mscheme_fact_GIIIc", &UnitTest::Mscheme_fact_GIIIc)
      .def("Mscheme_fact_GIVa", &UnitTest::Mscheme_fact_GIVa)
      .def("Mscheme_fact_GIVb_chi", &UnitTest::Mscheme_fact_GIVb_chi)
      .def("Mscheme_fact_GIVc", &UnitTest::Mscheme_fact_GIVc)
      .def("Mscheme_comm231tts_wick", &UnitTest::Mscheme_comm231tts_wick)
      .def("Mscheme_comm132tts_wick", &UnitTest::Mscheme_comm132tts_wick)
      .def("Mscheme_comm232tts_wick", &UnitTest::Mscheme_comm232tts_wick)
      .def("Mscheme_comm223tts_wick", &UnitTest::Mscheme_comm223tts_wick)


      ;

  //  py::class_<SymmMatrix<double>>(m,"SymmMatrix")
  //     .def(py::init<size_t>())
  //     .def(py::init<size_t,int>())
  //     .def("Get",&SymmMatrix<double>::Get)
  //     .def("Put",&SymmMatrix<double>::Put)
  //     .def("FullMatrix",&SymmMatrix<double>::FullMatrix)
  //  ;

  m.def("BuildVersion", version::BuildVersion);

  m.def("TCM_Op", imsrg_util::TCM_Op);
  m.def("Trel_Op", imsrg_util::Trel_Op);
  m.def("R2CM_Op", imsrg_util::R2CM_Op);
  m.def("HCM_Op", imsrg_util::HCM_Op);
  m.def("NumberOp", imsrg_util::NumberOp);
  m.def("RSquaredOp", imsrg_util::RSquaredOp);
  m.def("RpSpinOrbitCorrection", imsrg_util::RpSpinOrbitCorrection);
  m.def("E0Op", imsrg_util::E0Op);
  m.def("AllowedFermi_Op", imsrg_util::AllowedFermi_Op);
  m.def("AllowedGamowTeller_Op", imsrg_util::AllowedGamowTeller_Op);
  m.def("ElectricMultipoleOp", imsrg_util::ElectricMultipoleOp);
  m.def("MagneticMultipoleOp", imsrg_util::MagneticMultipoleOp);
  m.def("SchiffOp", imsrg_util::SchiffOp);
  m.def("Sigma_Op", imsrg_util::Sigma_Op);
  m.def("Isospin2_Op", imsrg_util::Isospin2_Op);
  m.def("LdotS_Op", imsrg_util::LdotS_Op);
  m.def("HO_density", imsrg_util::HO_density);
  m.def("GetOccupationsHF", imsrg_util::GetOccupationsHF);
  m.def("GetDensity", imsrg_util::GetDensity);
  m.def("Calculate_p1p2_all", imsrg_util::Calculate_p1p2_all);
  m.def("Single_Ref_1B_Density_Matrix",
        imsrg_util::Single_Ref_1B_Density_Matrix);
  m.def("Get_Charge_Density", imsrg_util::Get_Charge_Density);
  m.def("Embed1BodyIn2Body", imsrg_util::Embed1BodyIn2Body);
  m.def("RadialIntegral", imsrg_util::RadialIntegral);
  m.def("RadialIntegral_RpowK", imsrg_util::RadialIntegral_RpowK);
  m.def("RadialIntegral_Gauss", imsrg_util::RadialIntegral_Gauss, py::arg("na"),
        py::arg("la"), py::arg("nb"), py::arg("lb"), py::arg("sig"));
  m.def("RPA_resummed_1b", imsrg_util::RPA_resummed_1b, py::arg("OpIn"),
        py::arg("H"), py::arg("mode"));
  m.def("FirstOrderCorr_1b", imsrg_util::FirstOrderCorr_1b, py::arg("OpIn"),
        py::arg("H"));
  m.def("FrequencyConversionCoeff", imsrg_util::FrequencyConversionCoeff);
  m.def("OperatorFromString", imsrg_util::OperatorFromString);
  m.def("HO_Radial_psi", imsrg_util::HO_Radial_psi, py::arg("n"), py::arg("l"),
        py::arg("hw"), py::arg("r"));
  m.def("MBPT2_SpectroscopicFactor", imsrg_util::MBPT2_SpectroscopicFactor);
  m.def("SerberTypePotential", imsrg_util::SerberTypePotential,
        py::arg("modelspace"), py::arg("V0"), py::arg("mu"), py::arg("A"),
        py::arg("B"), py::arg("C"), py::arg("D"));

  m.def("CG", AngMom::CG);
  m.def("ThreeJ", AngMom::ThreeJ);
  m.def("SixJ", AngMom::SixJ);
  m.def("NineJ", AngMom::NineJ);
  m.def("NormNineJ", AngMom::NormNineJ);
  m.def("Moshinsky", AngMom::Moshinsky, py::arg("N"), py::arg("LAM"),
        py::arg("n"), py::arg("lam"), py::arg("n1"), py::arg("l1"),
        py::arg("n2"), py::arg("l2"), py::arg("L"), py::arg("BETA"));
  m.def("TalmiB", AngMom::TalmiB);
  m.def("TalmiI", imsrg_util::TalmiI);
  m.def("Tcoeff", AngMom::Tcoeff);
  m.def("SetUseGooseTank", BCH::SetUseGooseTank);
  m.def("SetUseIMSRG3", Commutator::SetUseIMSRG3);
  m.def("SetUseIMSRG3N7", Commutator::SetUseIMSRG3N7);
  m.def("FillFactorialLists", AngMom::FillFactorialLists);
  m.def("factorial", AngMom::factorial);
  m.def("double_fact", AngMom::double_fact);
  m.def("AngMomJmin", AngMom::Jmin);
  m.def("AngMomJmax", AngMom::Jmax);

  m.attr("HBARC") = py::float_(PhysConst::HBARC);
  m.attr("M_PROTON") = py::float_(PhysConst::M_PROTON);
  m.attr("M_NEUTRON") = py::float_(PhysConst::M_NEUTRON);
  m.attr("M_NUCLEON") = py::float_(PhysConst::M_NUCLEON);
  m.attr("M_ELECTRON") = py::float_(PhysConst::M_ELECTRON);
  m.attr("M_PION_CHARGED") = py::float_(PhysConst::M_PION_CHARGED);
  m.attr("M_PION_NEUTRAL") = py::float_(PhysConst::M_PION_NEUTRAL);
  m.attr("NUCLEON_VECTOR_G") = py::float_(PhysConst::NUCLEON_VECTOR_G);
  m.attr("NUCLEON_AXIAL_G") = py::float_(PhysConst::NUCLEON_AXIAL_G);
  m.attr("PROTON_SPIN_G") = py::float_(PhysConst::PROTON_SPIN_G);
  m.attr("NEUTRON_SPIN_G") = py::float_(PhysConst::NEUTRON_SPIN_G);
  m.attr("ELECTRON_SPIN_G") = py::float_(PhysConst::ELECTRON_SPIN_G);
  m.attr("ALPHA_FS") = py::float_(PhysConst::ALPHA_FS);
  m.attr("F_PI") = py::float_(PhysConst::F_PI);
  m.attr("HARTREE") = py::float_(PhysConst::HARTREE);
  m.attr("PROTON_RCH2") = py::float_(PhysConst::PROTON_RCH2);
  m.attr("NEUTRON_RCH2") = py::float_(PhysConst::NEUTRON_RCH2);
  m.attr("DARWIN_FOLDY") = py::float_(PhysConst::DARWIN_FOLDY);
}
