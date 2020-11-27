template <class T>
int setupAide::getArgs(string key, T& t){
  vector<T> m;

  getArgs(key,m);

  if(m.size()){
    t = m[0]; // TW

    return 1;
  }

  //printf("Failed to find [%s].\n", key.c_str());
  return 0;
}

template <class T>
int setupAide::getArgs(string key, vector<T>& m){
  stringstream args;
  vector<T> argv;
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
