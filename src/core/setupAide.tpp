template <class T>
int setupAide::getArgs(std::string key, T& t){
  std::vector<T> m;

  getArgs(key,m);

  if(m.size()){
    t = m[0]; // TW

    return 1;
  }

  //printf("Failed to find [%s].\n", key.c_str());
  return 0;
}

template <class T>
int setupAide::getArgs(std::string key, std::vector<T>& m){
  std::stringstream args;
  std::vector<T> argv;
  int argc;
  T input;

  args.str( getArgs(key) );

  while(args >> input)
    argv.push_back(input);

  argc = argv.size();

  if(!argc){
    //printf("Failed to find [%s].\n", key.c_str());
    return 0;
  }

  m.resize(argc);

  for(int i=0; i<argc; i++) // TW
    m[i] = argv[i];

  return 1;
}
