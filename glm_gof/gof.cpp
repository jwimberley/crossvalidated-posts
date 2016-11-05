Matrix InformationMatrix(const Table& table) { 
        
      int num = table[0].dpi.size();
      Matrix fisher(num,num);
        
      for (const auto& it : table) {
        double ppinv = it.pi*(1-it.pi);
        Matrix mtemp(num,num);
        for (unsigned int k = 0; k < num; k++) {
          for (unsigned int l = 0; l < num; l++) {
    	mtemp(k,l) = it.dpi[k]*it.dpi[l]/ppinv;
          }
        }
        fisher += it.W*mtemp;
      }
    
      return fisher;
    }
    
    Vector GenericScoreVector(const Table& table, std::function<double(bool,double)> metric) {
        
      int num = table[0].dpi.size();
      Vector score(num);
      
      for (const auto& it : table) {
        Vector grad(num);
        for (unsigned int k = 0; k < num; k++) {
          double pip = it.dpi[k];
          grad(k) = pip*(metric(true,it.pi)-metric(false,it.pi));
        }
        score += it.W*grad;
      }
    
      return score;
    }
    
    double GenericMetric(const Table& table, std::function<double(bool,double)> metric, bool selfCalibrated)
    {
      
      Matrix fisher = (selfCalibrated) ? InformationMatrix(table) : Matrix(1,1);
      Vector score = (selfCalibrated) ? GenericScoreVector(table, metric) : Vector(1);
      Vector scorep = (selfCalibrated) ? fisher.LinearSolve(score) : Vector(1); // fisher^-1 * score
      
      double difference = 0.0;
      double observed = 0.0;
      double expected = 0.0;
      double variance = 0.0; // to be reduced by a quadratic form
    
      for (const auto& it : table) 
        {
          double t = metric(true,it.pi);
          double f = metric(false,it.pi);
          double obs = (it.correct) ? t : f;
          double exp = it.pi*t + (1-it.pi)*f;
          double diff = obs - exp;
          double var = it.pi*(t-exp)*(t-exp) + (1-it.pi)*(f-exp)*(f-exp);
          if (selfCalibrated) {
            double a = t - f;
            double mu = scorep*it.dpi;
            var -= a*mu;
          }
          observed += it.W*obs;
          expected += it.W*exp;
          difference += it.W*diff;
          variance += it.W*var;
        }
      
      double z = (difference)/sqrt(std::abs(variance));
      return z;
    }
    
    
    double UngroupedDevianceTest(const Table& table, bool selfCalibrated) {
      
      auto deviance = [] (bool c, double p) {
        if (c)
          return -2*log(p);
        else
          return -2*log(1-p);
      };
      return GenericMetric(table,deviance,selfCalibrated);
    }
      
    double UngroupedPearsonTest(const Table& table, bool selfCalibrated) {
      auto pearson = [] (bool c, double p) {
        double den = p*(1-p);
        if (c)
          return (1-p)*(1-p)/den;
        else
          return p*p/den;
      };
      return GenericMetric(table,pearson,selfCalibrated);
    
    }
    
    double UngroupedSTest(const Table& table, bool selfCalibrated) {
      
      auto chch = [] (bool c, double p) {
        if (c)
          return (1-p)*(1-p);
        else
          return p*p;
      };
      return GenericMetric(table,chch,selfCalibrated);
    }
