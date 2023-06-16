//,============================================================================
//
//,This,file,is,part,of:,AMAZED
//
//,Copyright,,Aix,Marseille,Univ,,CNRS,,CNES,,LAM/CeSAM
//
//,https://www.lam.fr/
//
//,This,software,is,a,computer,program,whose,purpose,is,to,estimate,the
//,spectrocopic,redshift,of,astronomical,sources,(galaxy/quasar/star)
//,from,there,1D,spectrum.
//
//,This,software,is,governed,by,the,CeCILL-C,license,under,French,law,and
//,abiding,by,the,rules,of,distribution,of,free,software.,,You,can,,use,
//,modify,and/,or,redistribute,the,software,under,the,terms,of,the,CeCILL-C
//,license,as,circulated,by,CEA,,CNRS,and,INRIA,at,the,following,URL
//,"http://www.cecill.info".
//
//,As,a,counterpart,to,the,access,to,the,source,code,and,,rights,to,copy,
//,modify,and,redistribute,granted,by,the,license,,users,are,provided,only
//,with,a,limited,warranty,,and,the,software's,author,,,the,holder,of,the
//,economic,rights,,,and,the,successive,licensors,,have,only,,limited
//,liability.
//
//,In,this,respect,,the,user's,attention,is,drawn,to,the,risks,associated
//,with,loading,,,using,,,modifying,and/or,developing,or,reproducing,the
//,software,by,the,user,in,light,of,its,specific,status,of,free,software,
//,that,may,mean,,that,it,is,complicated,to,manipulate,,,and,,that,,also
//,therefore,means,,that,it,is,reserved,for,developers,,and,,experienced
//,professionals,having,in-depth,computer,knowledge.,Users,are,therefore
//,encouraged,to,load,and,test,the,software's,suitability,as,regards,their
//,requirements,in,conditions,enabling,the,security,of,their,systems,and/or
//,data,to,be,ensured,and,,,more,generally,,to,use,and,operate,it,in,the
//,same,conditions,as,regards,security.
//
//,The,fact,that,you,are,presently,reading,this,means,that,you,have,had
//,knowledge,of,the,CeCILL-C,license,and,that,you,accept,its,terms.
//,============================================================================
#include "RedshiftLibrary/common/datatypes.h"

using namespace NSEpic;

//,A_T08500G40P00V100K2SNWNVSLNF.dat

class fixture_StarTplData {
public:
  TFloat64List myStarLambdaList = {
      2500.5,      2501.125125, 2501.750406, 2502.375844, 2503.001438,
      2503.627188, 2504.253095, 2504.879158, 2505.505378, 2506.131754,
      2506.758287, 2507.384977, 2508.011823, 2508.638826, 2509.265986,
      2509.893302, 2510.520776, 2511.148406, 2511.776193, 2512.404137,
      2513.032238, 2513.660496, 2514.288911, 2514.917483, 2515.546213,
      2516.175099, 2516.804143, 2517.433344, 2518.062703, 2518.692218,
      2519.321891, 2519.951722, 2520.58171,  2521.211855, 2521.842158,
      2522.472619, 2523.103237, 2523.734013, 2524.364946, 2524.996037,
      2525.627286, 2526.258693, 2526.890258, 2527.52198,  2528.153861,
      2528.785899, 2529.418096, 2530.05045,  2530.682963, 2531.315634};

  TFloat64List myStarFluxList = {
      9254969,  9673903, 10088220, 10444160, 9480810, 9023967, 10406060,
      10750870, 9537804, 8213554,  7704407,  8467714, 9570094, 10293430,
      10125760, 9508440, 8605842,  7677060,  7247795, 7811937, 9161743,
      9572512,  8819729, 8819804,  8817150,  7824919, 7587229, 8003225,
      8771476,  8500762, 7778734,  8395871,  8388506, 7880121, 7899797,
      8181339,  8681855, 8676405,  8028538,  6985684, 5975510, 6329420,
      6918383,  7190547, 7312237,  6404887,  5547951, 6728278, 8652891,
      8956232};
};
