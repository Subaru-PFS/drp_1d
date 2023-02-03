// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/common/datatypes.h"

using namespace NSEpic;

static TFloat64List mySpectralList = {
    4.680282E+03, 4.680882E+03, 4.681482E+03, 4.682082E+03, 4.682682E+03,
    4.683282E+03, 4.683882E+03, 4.684482E+03, 4.685083E+03, 4.685683E+03,
    4.686283E+03, 4.686883E+03, 4.687483E+03, 4.688083E+03, 4.688683E+03,
    4.689283E+03, 4.689883E+03, 4.690483E+03, 4.691083E+03, 4.691683E+03,
    4.692283E+03, 4.692883E+03, 4.693483E+03, 4.694083E+03, 4.694684E+03,
    4.695283E+03, 4.695883E+03, 4.696483E+03, 4.697083E+03, 4.697684E+03,
    4.698284E+03, 4.698884E+03, 4.699484E+03, 4.700084E+03, 4.700684E+03,
    4.701284E+03, 4.701884E+03, 4.702484E+03, 4.703084E+03, 4.703684E+03,
    4.704284E+03, 4.704884E+03, 4.705484E+03, 4.706084E+03, 4.706685E+03,
    4.707285E+03, 4.707885E+03, 4.708485E+03, 4.709085E+03, 4.709685E+03,
    4.710285E+03, 4.710885E+03, 4.711485E+03, 4.712085E+03};

static TFloat64List myFluxList = {
    -1.059679E-18, 2.676246E-18, 5.185873E-18,  1.372732E-17,  4.928306E-18,
    2.443292E-19,  4.087630E-18, 3.933421E-19,  1.698866E-18,  1.999020E-18,
    3.494401E-18,  3.967686E-18, 5.091954E-18,  9.332690E-18,  5.365104E-18,
    7.323179E-18,  1.411062E-18, 3.868855E-18,  3.213580E-18,  1.532060E-18,
    2.229017E-18,  4.934030E-18, 5.108419E-18,  2.345380E-17,  6.194209E-17,
    5.760186E-17,  1.553542E-17, 2.149000E-18,  -3.494747E-18, 2.118534E-17,
    5.809586E-17,  4.255493E-17, 1.136803E-17,  1.756577E-18,  1.573434E-18,
    -3.550685E-18, 2.317096E-18, 2.476682E-18,  -2.245990E-18, 3.596576E-18,
    5.688406E-18,  3.259921E-18, 9.329072E-19,  -6.800624E-19, 1.774613E-18,
    4.460985E-18,  1.519357E-18, -2.147758E-18, -1.688685E-19, -5.715160E-19,
    8.158859E-19,  6.349788E-18, 4.712458E-18,  3.473506E-18};

static TFloat64List myNoiseList = {
    4.774209E-18, 4.805260E-18, 4.897541E-18, 5.022887E-18, 5.087175E-18,
    5.019457E-18, 4.869365E-18, 4.749990E-18, 4.693996E-18, 4.676343E-18,
    4.671902E-18, 4.669868E-18, 4.668044E-18, 4.666224E-18, 4.664510E-18,
    4.662815E-18, 4.661106E-18, 4.659397E-18, 4.657705E-18, 4.656002E-18,
    4.654295E-18, 4.652601E-18, 4.734505E-18, 5.274864E-18, 6.127447E-18,
    6.045540E-18, 4.901570E-18, 4.643016E-18, 4.740873E-18, 5.627134E-18,
    6.268595E-18, 5.567977E-18, 4.836006E-18, 4.632385E-18, 4.630944E-18,
    4.629538E-18, 4.628140E-18, 4.626769E-18, 4.625365E-18, 4.623963E-18,
    4.622573E-18, 4.621174E-18, 4.619774E-18, 4.618386E-18, 4.617014E-18,
    4.615645E-18, 4.614266E-18, 4.612899E-18, 4.611641E-18, 4.610969E-18,
    4.615064E-18, 4.637263E-18, 4.680976E-18, 4.718082E-18};

static TFloat64List mySpectralListLog = {
    4680,   4680.5, 4681,   4681.5, 4682,   4682.5, 4683,   4683.5, 4684,
    4684.5, 4685,   4685.5, 4686,   4686.5, 4687,   4687.5, 4688,   4688.5,
    4689,   4689.5, 4690,   4690.5, 4691,   4691.5, 4692,   4692.5, 4693,
    4693.5, 4694,   4694.5, 4695,   4695.5, 4696,   4696.5, 4697,   4697.5,
    4698,   4698.5, 4699,   4699.5, 4700,   4700.5, 4701,   4701.5, 4702,
    4702.5, 4703,   4703.5, 4704,   4704.5, 4705,   4705.5, 4706,   4706.5};
