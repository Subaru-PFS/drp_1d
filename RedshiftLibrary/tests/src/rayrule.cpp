#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/rule.h>
//#include <RedshiftLibrary/linemodel/multiline.h>
//#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/processflow/processflow.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/method/linemodelsolve.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

#include <math.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

using namespace std;
using namespace boost;

BOOST_AUTO_TEST_SUITE(RayRule)

class GoodDataConcreteRule : public CRule
{
public:
  Bool Check( CLineModelElementList& LineModelElementList );
  void SetUp( Bool EnabledArgument, ... );
private:
  void Correct( CLineModelElementList& LineModelElementList );
};

void GoodDataConcreteRule::SetUp( Bool EnabledArgument, ... )
{
}

Bool GoodDataConcreteRule::Check( CLineModelElementList& LineModelElementList )
{
  return true;
}

void GoodDataConcreteRule::Correct( CLineModelElementList& LineModelElementList )
{
  return;
}

class BadDataConcreteRule : public CRule
{
public:
  Bool Check( CLineModelElementList& LineModelElementList );
  void SetUp( Bool EnabledArgument, ... );
private:
  void Correct( CLineModelElementList& LineModelElementList );
};

void BadDataConcreteRule::SetUp( Bool EnabledArgument, ... )
{
}

Bool BadDataConcreteRule::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

void BadDataConcreteRule::Correct( CLineModelElementList& LineModelElementList )
{
  std::cout << "Correct" << std::endl;
  int dataSize = LineModelElementList.m_Elements.size( );
  std::cout << "Correct dataSize == " << dataSize << std::endl;
  for ( int i=0; i<dataSize; i++ )
    {
      std::cout << "Correcting LineModelElementList.m_Elements[" << i << "] from " << LineModelElementList.m_Elements[i]->GetFittedAmplitude(0);
      LineModelElementList.m_Elements[i]->SetFittedAmplitude ( -2, -2 );
      std::cout << " to " << LineModelElementList.m_Elements[i]->GetFittedAmplitude(0) << "." << std::endl;
    }
  return;
}

CLineModelElementList GetData ( void )
{
  std::string spc, noise;
  spc = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_1.fits";
  noise = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_1_noise.fits";
  CProcessFlowContext ctx;
  CProcessFlow processFlow;
  TFloat64Range redshiftRange = TFloat64Range( 0.0, 0.1 );
  TFloat64Range spcLambdaRange = TFloat64Range( 3800.0, 12000.0 );
  std::shared_ptr<CParameterStore> params = std::shared_ptr<CParameterStore>( new CParameterStore() );
  params->Set( "lambdaRange", spcLambdaRange);
  params->Set( "redshiftRange",  redshiftRange);
  params->Set( "redshiftStep", 0.01);
  params->Set( "smoothWidth", (Int64)0 );
  params->Set( "templateCategoryList", TStringList { "galaxy" } );
  params->Set( "method", "linemodel");

  std::string procID = "processing_id_unused";
  std::shared_ptr<CClassifierStore> classifStore = std::shared_ptr<CClassifierStore>( new CClassifierStore() );
  Bool retVal = ctx.Init( spc.c_str(), noise.c_str(), procID, NULL, "../test/data/LinemodelRulesTestCase/raycatalog_test_elratiorules.txt", params, classifStore );
  BOOST_CHECK( retVal == true );
  retVal = processFlow.Process( ctx ); // Segmentation fault
  BOOST_CHECK( retVal == true );
  // Create redshift initial list by spanning redshift acdross the given range, with the given delta
  Float64 redshiftStep = 0.01;
  TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
  Float64 resolution = 2350.0;
  Float64 emissionVelocity = 100.0;
  Float64 absorptionVelocity = 300.0;
  std::string opt_rules = "all";
  std::string opt_rigidity = "rules";
  std::string unused_calibrationPath="";


  //these tplcatalog related variables are unused here.
  CTemplateCatalog tplCatalog;
  Bool retValue = tplCatalog.Load( "../test/data/templatecatalog/" );
  TStringList tplCategories;

  //* Segmentation fault
  CLineModelElementList ReferenceArgument = CLineModelElementList ( ctx.GetSpectrum(),
								    ctx.GetSpectrumWithoutContinuum(),
                                                                    tplCatalog,
                                                                    tplCategories,
                                                                    unused_calibrationPath,
								    ctx.GetRayCatalog().GetList(),
								    std::string( "hybrid" ),
								    std::string( "fromspectrum" ),
								    std::string( "fixedvelocity" ),
								    resolution,
								    emissionVelocity,
								    absorptionVelocity,
                                                                    opt_rules,
                                                                    opt_rigidity);
  ReferenceArgument.LoadCatalog ( ctx.GetRayCatalog().GetList() );
  return ReferenceArgument;
}

Bool CompareData ( std::vector<boost::shared_ptr<CLineModelElement> > a,
		   std::vector<boost::shared_ptr<CLineModelElement> > b )
{
  int dataSize = a.size();
  std::cout << "CompareData: dataSize = " << dataSize << std::endl;
  for ( int i=0; i<dataSize; i++ )
    {
      auto aa = a [ i ]->GetFittedAmplitude( 0 );
      auto ab = b [ i ]->GetFittedAmplitude( 0 ); 
      std::cout << "CompareData: aa = " << aa << " ab = " << ab << endl;
      if ( aa != ab )
	{
	  return false;
	}
    }
  return true;
}

BOOST_AUTO_TEST_CASE(DisabledGoodCheck)
{
  // Disabled rule + Good data, Check() -> True.
  GoodDataConcreteRule aRule = GoodDataConcreteRule ( );
  CLineModelElementList someData = GetData ( );
  Bool test = aRule.Check ( someData );
  BOOST_CHECK( test == true );
}

BOOST_AUTO_TEST_CASE(DisabledBadCheck)
{
  // Disabled rule + Bad data, Check() -> False.
  BadDataConcreteRule aRule = BadDataConcreteRule ( );
  CLineModelElementList someData = GetData ( );
  Bool test = aRule.Check ( someData );
  BOOST_CHECK( test == false );
}

BOOST_AUTO_TEST_CASE(DisabledGoodApply)
{
  // Disabled rule + Good data, Apply() -> output == input.
  GoodDataConcreteRule aRule = GoodDataConcreteRule ( );
  CLineModelElementList someData = GetData ( );
  aRule.Apply ( someData );
  CLineModelElementList originalData = GetData ( );
  BOOST_CHECK( CompareData ( someData.m_Elements, originalData.m_Elements ) );
}

BOOST_AUTO_TEST_CASE(DisabledBadApply)
{
  // Disabled rule + Bad data, Apply() -> output == input.
  BadDataConcreteRule aRule = BadDataConcreteRule ( );
  CLineModelElementList someData = GetData ( );
  aRule.Apply ( someData );
  CLineModelElementList originalData = GetData ( );
  BOOST_CHECK( CompareData ( someData.m_Elements, originalData.m_Elements ) );
}

BOOST_AUTO_TEST_CASE(EnabledBadApply)
{
  // Enabled rule + Bad data, Apply() -> output, Check() -> True.
  BadDataConcreteRule aRule = BadDataConcreteRule ( );
  aRule.Enabled = true;
  CLineModelElementList someData = GetData ( );
  aRule.Apply ( someData );
  CLineModelElementList originalData = GetData ( );
  BOOST_CHECK( CompareData ( someData.m_Elements, originalData.m_Elements ) == false );
}

BOOST_AUTO_TEST_CASE(EnabledGoodApply)
{
  // Enabled rule + Good data, Apply() -> output == input.
  GoodDataConcreteRule aRule = GoodDataConcreteRule ( );
  aRule.Enabled = true;
  CLineModelElementList someData = GetData ( );
  aRule.Apply ( someData );
  CLineModelElementList originalData = GetData ( );
  BOOST_CHECK( CompareData ( someData.m_Elements, originalData.m_Elements ) );
}

BOOST_AUTO_TEST_SUITE_END()
