#include <time.h>
#include <iostream>
#include <istream>
#include <streambuf>
#include <sstream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>


#include<boost/tokenizer.hpp>
#include<string>

#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/ray/ray.h>



using namespace NSEpic;
struct membuf : std::streambuf
{
    membuf(char* begin, char* end) {
        this->setg(begin, begin, end);
    }
};

BOOST_AUTO_TEST_SUITE(test_line_model_result)
BOOST_AUTO_TEST_CASE(ResizeExtremaResults){
  CLineModelResult linemodelResult = CLineModelResult();

  linemodelResult.ResizeExtremaResults(2);
  BOOST_CHECK(linemodelResult.Extrema.size() == 2);
  BOOST_CHECK(linemodelResult.ExtremaMerit.size() == 2);
  BOOST_CHECK(linemodelResult.DeltaZ.size() == 2);
  BOOST_CHECK(linemodelResult.mTransposeM.size() == 2);
  BOOST_CHECK(linemodelResult.ExtremaLastPass.size() == 2);
  BOOST_CHECK(linemodelResult.Posterior.size() == 2);
  BOOST_CHECK(linemodelResult.StrongELSNR.size() == 2);
  BOOST_CHECK(linemodelResult.LogArea.size() == 2);
  BOOST_CHECK(linemodelResult.LogAreaCorrectedExtrema.size() == 2);
  BOOST_CHECK(linemodelResult.SigmaZ.size() == 2);
  BOOST_CHECK(linemodelResult.bic.size() == 2);
  BOOST_CHECK(linemodelResult.ContinuumIndexes.size() == 2);
  BOOST_CHECK(linemodelResult.OutsideLinesMask.size() == 2);
  BOOST_CHECK(linemodelResult.FittedTplName.size() == 2);
  BOOST_CHECK(linemodelResult.FittedTplAmplitude.size() == 2);
  BOOST_CHECK(linemodelResult.FittedTplDustCoeff.size() == 2);
  BOOST_CHECK(linemodelResult.FittedTplMeiksinIdx.size() == 2);

}

BOOST_AUTO_TEST_CASE(Save){
  CSpectrum s;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis( 15000, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<spectralAxis.GetSamplesCount(); k++){
    fluxAxis[k]=k;
  }

  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(15000);
  Float64* modelSamples = modelfluxAxis.GetSamples();
  for(Int32 k=0; k<modelfluxAxis.GetSamplesCount(); k++){
    modelSamples[k]=10.;
  }
  s.GetSpectralAxis() = spectralAxis;
  s.GetFluxAxis() = modelfluxAxis;
  CContinuumIndexes continuumIndexes;
  CContinuumIndexes::TContinuumIndexList indexesList = continuumIndexes.getIndexes( s, 0.0 );

  CMask mask = CMask();
  CLineModelResult linemodelResult = CLineModelResult();

  char buffer[] = "#\n0.1 \t 0.3 \t";

  membuf sbuf(buffer, buffer + sizeof(buffer));
  std::istream istream(&sbuf);
  linemodelResult.Load(istream);

  linemodelResult.ResizeExtremaResults(1);
  linemodelResult.Extrema[0]= 0.1;
  linemodelResult.ExtremaMerit[0]= 0.1;
  linemodelResult.DeltaZ[0]= 0.1;
  linemodelResult.mTransposeM[0]= 0.1;
  linemodelResult.ExtremaLastPass[0]= 0.1;
  linemodelResult.Posterior[0]= 0.1;
  linemodelResult.StrongELSNR[0]= 0.1;
  linemodelResult.LogArea[0]= 0.1;
  linemodelResult.LogAreaCorrectedExtrema[0]= 0.1;
  linemodelResult.SigmaZ[0]= 0.1;
  linemodelResult.bic[0]= 0.1;
  linemodelResult.ContinuumIndexes[0] = indexesList;
  linemodelResult.OutsideLinesMask[0]= mask;
  linemodelResult.FittedTplName[0]= "template";
  linemodelResult.FittedTplAmplitude[0]= 0.1;
  linemodelResult.FittedTplDustCoeff[0]= 0.1;
  linemodelResult.FittedTplMeiksinIdx[0]= 1;
  linemodelResult.dTransposeDNocontinuum =0.1;
  linemodelResult.dTransposeD = 0.1;

  NSEpic::COperatorResultStore resultStore = NSEpic::COperatorResultStore();
  NSEpic::CParameterStore parameStore =NSEpic::CParameterStore();
  CDataStore store(resultStore, parameStore);
  //
  std::ostringstream stream;
  // std::stringbuf str;
  // stream.rdbuf(&str);
  linemodelResult.Save(store, stream);
  std::ostringstream result;
  result << "#Redshifts\tChiSquare";
  result<< std::endl <<"0.1\t2.99999999999999988897769753748435e-01";
  result<< std::endl <<"#Extrema for z = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#ExtremaMerit for z = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#ExtremaDeltaZ for z = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#mTransposeM for z = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#ExtremaLastPass for z = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#BIC for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#POSTERIOR for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#SigmaZ for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#LogArea for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#ContinuumIndexes Color for each extrema = {<-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t>}";
  result<< std::endl <<"#ContinuumIndexes Break for each extrema = {<-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t-0.00000000000000000000000000000000\t>}";
  result<< std::endl <<"#StrongELSNR for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#dTransposeDNocontinuum for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#dTransposeD for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#FittedTplName for each extrema = {template\t}";
  result<< std::endl <<"#FittedTplAmplitude for each extrema = {0.10000000000000000555111512312578\t}";
  result<< std::endl <<"#FittedTplDustCoeff for each extrema = {0.100\t}";
  result<< std::endl <<"#FittedTplMeiksinIdx for each extrema = {1\t}";
  result<< std::endl <<"#FittedTplcorrTplName for each extrema = {\t}";
  result<< std::endl;

  // std::stringstream ss(result.str());
  // std::stringstream ss2(stream.str());
  // std::string to;
  // std::string to2;
  // while(std::getline(ss,to,'\n')){
  //   std::getline(ss2,to2,'\n');
  //   if(to == to2){
  //     BOOST_TEST_MESSAGE(""<<to <<";true");
  //   }else{
  //     BOOST_TEST_MESSAGE(""<<to <<";false");
  //   }
  // }
  BOOST_CHECK(stream.str() == result.str());
  // BOOST_TEST_MESSAGE("Save string:"<<stream.str()<<";");
  // BOOST_TEST_MESSAGE("result string:"<<result.str())<<";";
}

BOOST_AUTO_TEST_CASE(SaveLine){

  CLineModelResult linemodelResult = CLineModelResult();
  linemodelResult.Redshifts.resize(1);
  linemodelResult.Redshifts[0] = 0.1;

  NSEpic::COperatorResultStore resultStore = NSEpic::COperatorResultStore();
  NSEpic::CParameterStore parameStore =NSEpic::CParameterStore();
  CDataStore store(resultStore, parameStore);
  //
  std::ostringstream stream;
  linemodelResult.SaveLine(store, stream);
  std::ostringstream result;
  result << "LineModelResult\t1"<< std::endl;
  BOOST_CHECK(stream.str() == result.str());
}

BOOST_AUTO_TEST_CASE(GetNLinesOverCutThreshold){
  CRay ray = CRay("Abs",5500, 1, "SYM", 2, 10.2, 10.3, 10.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CRay ray2 = CRay("Abs",5500, 2, "SYM", 2, 10.2, 10.3, 10.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CRay ray3 = CRay("Abs",5500, 2, "SYM", 1, 10.2, 10.3, 10.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  //ray, ray2, ray3 are the same but not strong or emision
  CRay ray4 = CRay("Em",5520, 2, "SYM", 2, 10.2, 10.3, 20.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CLineModelResult linemodelResult = CLineModelResult();
  linemodelResult.Redshifts.push_back(0.6);
  linemodelResult.Redshifts.push_back(0.8);
  linemodelResult.Extrema.push_back(0.8);
  linemodelResult.Extrema.push_back(0.6);

  CLineModelResult::SLineModelSolution s1;
  s1.ElementId.push_back(0);
  s1.ElementId.push_back(1);
  s1.ElementId.push_back(2);
  s1.ElementId.push_back(1);
  s1.Amplitudes.push_back(2.1);
  s1.Amplitudes.push_back(2.1);
  s1.Amplitudes.push_back(2.1);
  s1.Amplitudes.push_back(1.2);
  s1.Rays.push_back(ray);
  s1.Rays.push_back(ray2);
  s1.Rays.push_back(ray3);
  s1.Rays.push_back(ray4);
  s1.Errors.push_back(0.5);
  s1.Errors.push_back(0.5);
  s1.Errors.push_back(0.5);
  s1.Errors.push_back(0.3);
  s1.FittingError.push_back(1.5);
  s1.FittingError.push_back(1.5);
  s1.FittingError.push_back(1.5);
  s1.FittingError.push_back(2.5);
  linemodelResult.LineModelSolutions.push_back(s1);

  CLineModelResult::SLineModelSolution s2;
  s2.ElementId.push_back(0);
  s2.ElementId.push_back(1);
  s2.ElementId.push_back(2);
  s2.ElementId.push_back(3);
  s2.Amplitudes.push_back(3.1);
  s2.Amplitudes.push_back(3.1);
  s2.Amplitudes.push_back(3.1);
  s2.Amplitudes.push_back(1.2);
  s2.Rays.push_back(ray);
  s2.Rays.push_back(ray2);
  s2.Rays.push_back(ray3);
  s2.Rays.push_back(ray4);
  s2.Errors.push_back(0.5);
  s2.Errors.push_back(0.5);
  s2.Errors.push_back(0.5);
  s2.Errors.push_back(0.3);
  s2.FittingError.push_back(1.5);
  s2.FittingError.push_back(1.5);
  s2.FittingError.push_back(1.5);
  s2.FittingError.push_back(2.5);
  linemodelResult.LineModelSolutions.push_back(s2);

//   A	E	FE	snr	fsnr
// S1.1	2,1	0,5	1,5	4,2	1,4
// S1.2	1,2	0,3	2,5	4	0,48
// S2.1	3,1	0,5	1,5	6,2	2,06666666666667
// S2.1	1,2	0,3	2,5	4	0,48


  BOOST_CHECK_EQUAL(1, linemodelResult.GetNLinesOverCutThreshold(0,5., 2.));
  BOOST_CHECK_EQUAL(0, linemodelResult.GetNLinesOverCutThreshold(0,5., 2.1));
  BOOST_CHECK_EQUAL(1, linemodelResult.GetNLinesOverCutThreshold(1,4.1, 1.3));
  BOOST_CHECK_EQUAL(0, linemodelResult.GetNLinesOverCutThreshold(1,4.1, 1.5));//fsnr to high
  BOOST_CHECK_EQUAL(1, linemodelResult.GetNLinesOverCutThreshold(1,3.9, 0.4));//ray2 and ray 4 same element ID

}

BOOST_AUTO_TEST_CASE(GetExtremaMerit){
  CLineModelResult linemodelResult = CLineModelResult();
  linemodelResult.ResizeExtremaResults(2);
  linemodelResult.Extrema[0]= 0.1;
  linemodelResult.Extrema[1]= 0.3;
  linemodelResult.ExtremaMerit[0]= 0.2;
  linemodelResult.ExtremaMerit[1]= 0.4;

  BOOST_CHECK_CLOSE(linemodelResult.GetExtremaMerit(1), 0.4,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.GetExtremaMerit(0), 0.2,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.GetExtremaMerit(2), -1,1e-12);

}

BOOST_AUTO_TEST_CASE(SaveFloatList){
  CLineModelResult linemodelResult = CLineModelResult();
  TFloat64List fl1;
  fl1.push_back(0.2);
  fl1.push_back(0.3);

  TFloat64List fl2;
  std::ostringstream result;
  result <<"blabla 0.2\t0.3\t}"<<std::endl;

  std::ostringstream result2;
  result2 <<"blabla }"<<std::endl;



  std::ostringstream stream1;
  linemodelResult.SaveFloatList(stream1, fl1, "blabla ", true);
  BOOST_CHECK_EQUAL(stream1.str(),result.str());

  std::ostringstream stream2;
  linemodelResult.SaveFloatList(stream2, fl1, "blabla ", false);
  BOOST_CHECK_EQUAL(stream2.str(),result.str());

  std::ostringstream stream3;
  linemodelResult.SaveFloatList(stream3, fl2, "blabla ", false);
  BOOST_CHECK_EQUAL(stream3.str(),result2.str());

  std::ostringstream stream4;
  linemodelResult.SaveFloatList(stream4, fl2, "blabla ", true);
  BOOST_CHECK_EQUAL(stream4.str(),"");
}
BOOST_AUTO_TEST_CASE(Load){


  char buffer[] = "#\n0.1 \t 0.3 \t\n0.2 \t 1 \t";

  membuf sbuf(buffer, buffer + sizeof(buffer));
  std::istream stream(&sbuf);

  CLineModelResult linemodelResult = CLineModelResult();
  linemodelResult.Load(stream);
  BOOST_CHECK_EQUAL(2, linemodelResult.Redshifts.size());
  BOOST_CHECK_CLOSE(linemodelResult.Redshifts[0], 0.1,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.Redshifts[1], 0.2,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.ChiSquare[0], 0.3,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.ChiSquare[1], 1,1e-12);

  char buffer2[] ="#blabla\n0.1 \t 0.3 \t\n0.2 \t 0.4 \t";

  membuf sbuf2(buffer2, buffer2 + sizeof(buffer2));
  std::istream stream2(&sbuf2);
  linemodelResult.Load(stream2);
  BOOST_CHECK_EQUAL(0, linemodelResult.Redshifts.size());


  char buffer3[] ="#\n0.5 \t 0.7 \t\nblabla0.2 \t 0.4 \t";

  membuf sbuf3(buffer3, buffer3 + sizeof(buffer3));
  std::istream stream3(&sbuf3);
  linemodelResult.Load(stream3);
  BOOST_CHECK_EQUAL(1, linemodelResult.Redshifts.size());
  BOOST_CHECK_CLOSE(linemodelResult.Redshifts[0], 0.5,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.ChiSquare[0], 0.7,1e-12);

  char buffer4[] ="#\n1.1 \t 1.3 \t\n0.2 \t blabla0.4 \t";

  membuf sbuf4(buffer4, buffer4 + sizeof(buffer4));
  std::istream stream4(&sbuf4);
  linemodelResult.Load(stream4);
  BOOST_CHECK_EQUAL(1, linemodelResult.Redshifts.size());
  BOOST_CHECK_CLOSE(linemodelResult.Redshifts[0], 1.1,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.ChiSquare[0], 1.3,1e-12);


  char buffer5[] ="#\n2.1 \t 2.3 \t\n0.2 ";

  membuf sbuf5(buffer5, buffer5 + sizeof(buffer5));
  std::istream stream5(&sbuf5);
  linemodelResult.Load(stream5);
  BOOST_CHECK_EQUAL(1, linemodelResult.Redshifts.size());
  BOOST_CHECK_CLOSE(linemodelResult.Redshifts[0], 2.1,1e-12);
  BOOST_CHECK_CLOSE(linemodelResult.ChiSquare[0], 2.3,1e-12);
}
BOOST_AUTO_TEST_SUITE_END()
