void printHeader()
{
  cout << R"(                 __    ____  _____)" << endl
       << R"(   ____   ___   / /__ / __ \/ ___/)" << endl
       << R"(  / __ \ / _ \ / //_// /_/ /\__ \ )" << endl
       << R"( / / / //  __// ,<  / _, _/___/ / )" << endl
       << R"(/_/ /_/ \___//_/|_|/_/ |_|/____/  )"
       << "v" << NEKRS_VERSION << "." << NEKRS_SUBVERSION << "." << NEKRS_PATCHVERSION << GITCOMMITHASH << endl
       << endl
       << "COPYRIGHT (c) 2019-2021 UCHICAGO ARGONNE, LLC" << endl
       << endl;
}
