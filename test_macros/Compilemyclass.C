void Compilemyclass(){
  gROOT -> ProcessLineSync(".L ../MathFuncsLib.cxx+") ;
  gROOT -> ProcessLineSync(".x ../Binning.cxx+") ;
}
