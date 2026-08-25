void parseRadiationSection(const int rank, setupAide &options, inipp::Ini *ini)
{
  if (!ini->sections.count("radiation")) {
    options.setArgs("RADIATION", "FALSE");
    return;
  }
  options.setArgs("RADIATION", "TRUE");

  std::string sbuf;

  if (ini->extract("radiation", "radiatingboundaryids", sbuf)) {
    options.setArgs("RADIATION RADIATING BOUNDARY IDS", sbuf);
  } else {
    append_error("[RADIATION] requires radiatingBoundaryIDs\n");
  }

  if (ini->extract("radiation", "obstructionboundaryids", sbuf)) {
    options.setArgs("RADIATION OBSTRUCTION BOUNDARY IDS", sbuf);
  }

  int nSamples = 4096;
  ini->extract("radiation", "nsamples", nSamples);
  options.setArgs("RADIATION NSAMPLES", std::to_string(nSamples));

  int seed = 0;
  ini->extract("radiation", "seed", seed);
  options.setArgs("RADIATION SEED", std::to_string(seed));

  std::string writeMatrix = "true";
  ini->extract("radiation", "writematrix", writeMatrix);
  options.setArgs("RADIATION WRITE MATRIX", checkForTrue(lowerCase(writeMatrix)) ? "TRUE" : "FALSE");

  std::string outputFile;
  if (ini->extract("radiation", "outputfile", outputFile)) {
    options.setArgs("RADIATION OUTPUT FILE", outputFile);
  }

  std::string cache = "true";
  ini->extract("radiation", "cache", cache);
  options.setArgs("RADIATION CACHE", checkForTrue(lowerCase(cache)) ? "TRUE" : "FALSE");

  if (ini->extract("radiation", "emissivity", sbuf)) {
    options.setArgs("RADIATION EMISSIVITY", sbuf);
  }

  double stefanBoltzmann = 5.670374419e-8;
  ini->extract("radiation", "stefanboltzmann", stefanBoltzmann);
  options.setArgs("RADIATION STEFAN BOLTZMANN", to_string_f(stefanBoltzmann));

  int updateFrequency = 10;
  ini->extract("radiation", "updatefrequency", updateFrequency);
  options.setArgs("RADIATION UPDATE FREQUENCY", std::to_string(updateFrequency));

  double radiosityTolerance = 1e-6;
  ini->extract("radiation", "radiositytolerance", radiosityTolerance);
  options.setArgs("RADIATION RADIOSITY TOLERANCE", to_string_f(radiosityTolerance));

  int radiosityMaxIters = 100;
  ini->extract("radiation", "radiositymaxiters", radiosityMaxIters);
  options.setArgs("RADIATION RADIOSITY MAX ITERS", std::to_string(radiosityMaxIters));
}
