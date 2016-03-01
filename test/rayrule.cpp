#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/rule.h>
//#include <epic/redshift/linemodel/multiline.h>
//#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/processflow/processflow.h>
#include <epic/redshift/processflow/parameterstore.h>
#include <epic/redshift/processflow/context.h>
#include <epic/redshift/method/linemodelsolve.h>
#include <epic/redshift/method/linemodelsolveresult.h>
#include <epic/redshift/linemodel/modelfittingresult.h>

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
  Bool Check( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements );
  void SetUp( Bool EnabledArgument, ... );
private:
  void Correct( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements );
};

void GoodDataConcreteRule::SetUp( Bool EnabledArgument, ... )
{
}

Bool GoodDataConcreteRule::Check( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  return true;
}

void GoodDataConcreteRule::Correct( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  return;
}

class BadDataConcreteRule : public CRule
{
public:
  Bool Check( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements );
  void SetUp( Bool EnabledArgument, ... );
private:
  void Correct( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements );
};

void BadDataConcreteRule::SetUp( Bool EnabledArgument, ... )
{
}

Bool BadDataConcreteRule::Check( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  return false;
}

void BadDataConcreteRule::Correct( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  std::cout << "Correct" << std::endl;
  int dataSize = LinemodelElements.size( );
  std::cout << "Correct dataSize == " << dataSize << std::endl;
  for ( int i=0; i<dataSize; i++ )
    {
      std::cout << "Correcting LinemodelElements[" << i << "] from " << LinemodelElements[i]->GetFittedAmplitude(0);
      LinemodelElements[i]->SetFittedAmplitude ( -2, -2 );
      std::cout << " to " << LinemodelElements[i]->GetFittedAmplitude(0) << "." << std::endl;
    }
  return;
}

void GetData ( std::vector<boost::shared_ptr<CLineModelElement> >& ReferenceArgument )
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
  Bool retVal = ctx.Init( spc.c_str(), noise.c_str(), NULL, "../test/data/LinemodelRulesTestCase/raycatalog_test_elratiorules.txt",params );
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
  //* Segmentation fault
  CLineModelElementList lineModel = CLineModelElementList ( ctx.GetSpectrum(),
							    ctx.GetSpectrumWithoutContinuum(),
							    ctx.GetRayCatalog().GetList(),
							    std::string( "hybrid" ),
							    std::string( "fromspectrum" ),
							    std::string( "fixedvelocity" ),
							    resolution,
							    emissionVelocity,
							    absorptionVelocity,
							    opt_rules );
  lineModel.LoadCatalog ( ctx.GetRayCatalog().GetList() );
  ReferenceArgument = lineModel.m_Elements;
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
  std::vector<boost::shared_ptr<CLineModelElement> > someData;
  GetData ( someData );
  Bool test = aRule.Check ( someData );
  BOOST_CHECK( test == true );
}

BOOST_AUTO_TEST_CASE(DisabledBadCheck)
{
  // Disabled rule + Bad data, Check() -> False.
  BadDataConcreteRule aRule = BadDataConcreteRule ( );
  std::vector<boost::shared_ptr<CLineModelElement> > someData;
  GetData ( someData );
  Bool test = aRule.Check ( someData );
  BOOST_CHECK( test == false );
}

BOOST_AUTO_TEST_CASE(DisabledGoodApply)
{
  // Disabled rule + Good data, Apply() -> output == input.
  GoodDataConcreteRule aRule = GoodDataConcreteRule ( );
  std::vector<boost::shared_ptr<CLineModelElement> > someData;
  GetData ( someData );
  aRule.Apply ( someData );
  std::vector<boost::shared_ptr<CLineModelElement> > originalData;
  GetData ( originalData );
  BOOST_CHECK( CompareData ( someData, originalData ) );
}

BOOST_AUTO_TEST_CASE(DisabledBadApply)
{
  // Disabled rule + Bad data, Apply() -> output == input.
  BadDataConcreteRule aRule = BadDataConcreteRule ( );
  std::vector<boost::shared_ptr<CLineModelElement> > someData;
  GetData ( someData );
  aRule.Apply ( someData );
  std::vector<boost::shared_ptr<CLineModelElement> > originalData;
  GetData ( originalData );
  BOOST_CHECK( CompareData ( someData, originalData ) );
}

BOOST_AUTO_TEST_CASE(EnabledBadApply)
{
  // Enabled rule + Bad data, Apply() -> output, Check() -> True.
  BadDataConcreteRule aRule = BadDataConcreteRule ( );
  aRule.Enabled = true;
  std::vector<boost::shared_ptr<CLineModelElement> > someData;
  GetData ( someData );
  aRule.Apply ( someData );
  std::vector<boost::shared_ptr<CLineModelElement> > originalData;
  GetData ( originalData );
  BOOST_CHECK( CompareData ( someData, originalData ) == false );
}

BOOST_AUTO_TEST_CASE(EnabledGoodApply)
{
  // Enabled rule + Good data, Apply() -> output == input.
  GoodDataConcreteRule aRule = GoodDataConcreteRule ( );
  aRule.Enabled = true;
  std::vector<boost::shared_ptr<CLineModelElement> > someData;
  GetData ( someData );
  aRule.Apply ( someData );
  std::vector<boost::shared_ptr<CLineModelElement> > originalData;
  GetData ( originalData );
  BOOST_CHECK( CompareData ( someData, originalData ) );
}

BOOST_AUTO_TEST_SUITE_END()
