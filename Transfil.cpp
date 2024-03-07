#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


class Data {
public:
  int time = 0;
  
  // Host
  int pop_size;
  double tau = 0.00167; // death rate
  NumericVector bite_risk; // individual bite risk
  
  IntegerVector bednet; // 1 if using bednet, otherwise 0
  double bednet_coverage = 0.0;
  double bednet_reduction = 0.97;
  
  double k; 
  double min_k = 0.01;
  double max_k = 1.6;
  
  double importation_rate;
  double min_importation_rate = 0;
  double max_importation_rate = 0.0025;
  //double max_importation_rate = 0.005;
  
  // Vector
  double lambda = 10.0; // bite rate mosq/month
  // proportion of mosquitos which pick up infection when biting an infected host
  double g = 0.37; 
  double larvae = 5.0; // average l3 density
  double sigma = 5.0; // mosquito death rate
  //double kappa = 4.395; // L3 update and development for anophales
  double kappa; // L3 update and development for anophales
  double min_kappa = 4.395;
  double max_kappa = 4.395;
  double r; // L3 update and development for anophales
  double min_r = 0.055;
  double max_r = 0.055;
  int species = 1; // 0 for culex, 1 for anophales
  double s_N = 0.03; // probability of successful feeding
  double d_N = 0.41; // probability of death
  
  // Worm
  double psi_1 = 0.414; // Proportion L3 leaving mosquito per bite
  double psi_2 = 0.32; // Proportion of L3 leaving mosquito that enter host
  double s_2 = 0.00275; // Proportion of L3 entering host that developintoadults
  double mu = 0.0104; // worm death rate
  double gamma = 0.1; // mf death rate
  double alpha; // mf birth rate per female worm
  double min_alpha = 0.58;
  double max_alpha = 0.58;
  
  double vth; // vector to host ratio
  double min_vth = 1.0;
  double max_vth = 180.0;
  
  // individual level worm/mf burden
  IntegerVector WF; // female worms
  IntegerVector WM; // male worms
  NumericVector Mf; // individual mf density / 20 ml blood
  NumericVector age;
  
  // MDA
  NumericVector treatment_prob; // Probability of receiving treatment
  IntegerVector ttreat; // Time since treatment
  IntegerVector excluded;
  double treatment_coverage = 0.01; // Treatment coverage
  double systematic_compliance = 0.35; // Correlation between rounds
  int fecundity = 9; // month of worm infertility due to MDA
  double sys_ex = 0.05; // Proportion of the population that are excluded
  double prop_mf_killed = 0.99;
  double prop_worm_killed = 0.35;
  
  IntegerVector MDA_time; // Times to apply MDA
  NumericVector MDA_coverage; // Coverage of MDA
  
  // bednet settings
  IntegerVector bednet_time; // times to adjust bednets
  NumericVector bednet_coverage_vec; // bednet coverage
  
  // prevalence vectors
  NumericVector worm_prev;
  NumericVector mf_prev;
  NumericVector worm_prev67;
  NumericVector mf_prev67;
  IntegerVector pop67;
  
  // Constructor
  Data(double n_host_);
  Data(){};
  

  
  // Methods
  void timestep(double n_months);
  void burn_in();
  void worm_dynamics();
  void mf_dynamics();
  void larvae_dynamics();
  void calculate_prevalence();
  void aging();
  void importation();
  void set_MDA(List MDA_data);
  void set_bednet(List bednet_data);
  void run_MDA();
  void update_bednet();
  void set_parameter(std::string parameter, double value);
  List save();
  void load(List data);
};

// Constructor
Data::Data(double pop_size_) {
  pop_size = int(pop_size_);
  
  IntegerVector WF_(pop_size); 
  IntegerVector WM_(pop_size); 
  NumericVector Mf_(pop_size); 
  WF = WF_;
  WM = WM_;
  Mf = Mf_;
  
  // Simulate a random age between 0 and 100 years
  age = runif(pop_size, 0, 100*12);
  age = floor(age);
  
  // Draw the bite risk parameter
  k = R::runif(min_k, max_k);
  
  // Draw probabilities of treatment
  double beta_shape1 = 
    treatment_coverage*(1-systematic_compliance)/systematic_compliance;
  double beta_shape2 = 
    (1-treatment_coverage)*(1-systematic_compliance)/systematic_compliance;
  treatment_prob = rbeta(pop_size, beta_shape1, beta_shape2);
  
  excluded = rbinom(pop_size, 1, sys_ex);
  IntegerVector ttreat_(pop_size, 15*12);
  ttreat = ttreat_;
  bednet = rbinom(pop_size, 1, bednet_coverage);
  bite_risk = rgamma(pop_size, k, 1/k);
  vth = R::runif(min_vth, max_vth);
  r = R::runif(min_r,max_r);
  importation_rate = R::runif(min_importation_rate, max_importation_rate);
  alpha = R::runif(min_alpha, max_alpha);
  kappa = R::runif(min_kappa, max_kappa);
}



// expose the class to R
RCPP_MODULE(lf_module) {
  class_<Data>("Data")
  
  .constructor<double>()
  .constructor()
  
  .field("pop_size", &Data::pop_size)
  .field("worm_prev", &Data::worm_prev)
  .field("worm_prev67", &Data::worm_prev67)
  .field("mf_prev", &Data::mf_prev)
  .field("pop67", &Data::pop67)
  .field("vth", &Data::vth)
  .field("k", &Data::k)
  .field("r", &Data::r)
  .field("species", &Data::species)
  .field("Mf", &Data::Mf)
  .field("importation_rate", &Data::importation_rate)
  .field("alpha", &Data::alpha)
  .field("kappa", &Data::kappa)
  .field("time", &Data::time)
  
  .method("timestep", &Data::timestep)
  .method("burn_in", &Data::burn_in)
  .method("set_MDA", &Data::set_MDA)
  .method("set_bednet", &Data::set_bednet)
  .method("set_parameter", &Data::set_parameter)
  .method("save", &Data::save)
  .method("load", &Data::load)
  ;
}

void Data::timestep(double n_months) {
  //if(worm_prev.length() == 0) n_months++;
  for(int i = 0; i < n_months; i++) {
    run_MDA();
    update_bednet();
    worm_dynamics();
    mf_dynamics();
    larvae_dynamics();
    calculate_prevalence();
    aging();
    importation();
  }
}


void Data::worm_dynamics() {
  bool verbose = false;
  if(verbose) Rcout << "Updating worm dynamics\n";
  for(int i = 0; i < pop_size; i++) {
    double age_prop = age.at(i)/(9.0*12.0);
    double bite_rate = std::min(age_prop, 1.0);
    double BNR = bednet.at(i) ? s_N : 1.0;
    double birth_rate = 0.5 * lambda * bite_risk.at(i) * vth * psi_1 * psi_2 * 
      s_2 * bite_rate * larvae * BNR;
    if(verbose) {
      Rcout << "i = " << i << ", birth_rate = " << birth_rate << "\n";
      Rcout << " lambda = " << lambda << ", bite_risk = " << bite_risk.at(i)
            << ", vth = " << vth << ", psi_1 = " << psi_1
            << ", psi_2 = " << psi_2 << ", s_2 = " << s_2
            << ", bite_rate = " << bite_rate << ", larvae = " << larvae
            << ", BNR = " << BNR << "\n";
      Rcout << "Mean bite risk = " << mean(bite_risk) << "\n";
     }
    
    double male_death_rate = WM.at(i)*mu;
    double female_death_rate = WF.at(i)*mu;

    int male_births = (birth_rate > 0) ? R::rpois(birth_rate) : 0;
    int female_births = (birth_rate > 0) ? R::rpois(birth_rate) : 0;
    int male_deaths = (male_death_rate > 0) ? R::rpois(male_death_rate) : 0;
    int female_deaths = 
      (female_death_rate > 0) ? R::rpois(female_death_rate) : 0;
    
    WM.at(i) += male_births - male_deaths;
    WF.at(i) += female_births - female_deaths;
    if(WM.at(i) < 0) WM.at(i) = 0;
    if(WF.at(i) < 0) WF.at(i) = 0;
  }
}

void Data::mf_dynamics() {
  bool verbose = false;
  if(verbose) {
    Rcout << "Updating MF dynamics with alpha = " << alpha 
          << " and gamma = " << gamma << "\n";
  }
  for(int i = 0; i < pop_size; i++) {
    Mf.at(i) *= 1-gamma; 
    if(WM.at(i) == 0) continue;
    if(ttreat.at(i) <= fecundity) continue;
    Mf.at(i) += alpha * WF.at(i);
    if(Mf.at(i) < 0) Mf.at(i) = 0;
  }
}

void Data::larvae_dynamics() {
  bool verbose = false;
  NumericVector L3;
  switch(species) {
  case 0:
    //Culex
    L3 = kappa*(1-exp(-Mf*r/kappa));
  case 1:
    //Anophales
    L3 = kappa*(1-exp(-Mf*r/kappa))*(1-exp(-Mf*r/kappa));
  }
  double mean_larvae = sum(L3*bite_risk)/sum(bite_risk);
  double lambda_bednet = lambda*(bednet_coverage*(s_N-1) + 1);
  double sigma_bednet = sigma + sigma*d_N*bednet_coverage;
  larvae = lambda_bednet*g*mean_larvae/(sigma_bednet + lambda_bednet*psi_1);
  if(verbose) {
    Rcout << "Updating larvae dynamics with lambda_bednet = " << lambda_bednet
          << ", g = "<< g << ", mean_larvae = " << mean_larvae
          << ", sigma_bednet = " << sigma_bednet << ", psi_1 = " << psi_1
          << " - lambda = " << lambda << ", bednet_coverage = "
          << bednet_coverage << ", s_N = " << s_N
          << "\n";
  }
}

void Data::calculate_prevalence() {
  bool verbose = false;
  double worm_pos = 0;
  double worm67_pos = 0;
  double mf67_pos = 0;
  double cur_pop67 = 0;
  for(int i = 0; i < pop_size; i++) {
    double age_i = age.at(i);
    if(age_i < 12*5) continue;
    if(WM.at(i) > 0 || WF.at(i) > 0) worm_pos++;
    
    // Isolate ages 6-7 for TAS
    if(age_i < 12*6) continue;
    if(age_i > 12*8) continue;
    cur_pop67++;
    if(WM.at(i) > 0 || WF.at(i) > 0) worm67_pos++;
    if(Mf.at(i) > 0) mf67_pos++;
  }
  
  LogicalVector mf_ind = runif(pop_size, 0.0, 1.0) < 1-exp(-1*Mf);
  
  double current_wormprev67 = worm67_pos == 0 ? 0 : worm67_pos/cur_pop67;
  double current_mfprev67 = mf67_pos == 0 ? 0 : mf67_pos/cur_pop67;
  double current_wormprev = worm_pos == 0 ? 0 : worm_pos/pop_size;
  worm_prev.push_back(current_wormprev);
  mf_prev.push_back(mean(mf_ind));
  worm_prev67.push_back(current_wormprev67);
  mf_prev67.push_back(current_mfprev67);
  pop67.push_back(cur_pop67);
  

  // adjust importation rate if MDA has happened
  if(!is_true(any(MDA_time == time))) return;
  int prev_len = worm_prev.length();
  double prev_before = worm_prev[prev_len-2];
  if(prev_before == 0) return;
  double prevalence_ratio = worm_prev[prev_len-1]/prev_before;
  if(verbose) {
    Rcout << "Prevalence ratio = " << prevalence_ratio
          << ", prev before = " << worm_prev[prev_len-2]
          << ", prev after = " << worm_prev[prev_len-1]
          << ", importation rate before = " << importation_rate
          << ", importation rate after = " << importation_rate*prevalence_ratio
          << "\n";
  }
  importation_rate *= prevalence_ratio;
}

void Data::aging() {
  age = age + 1;
  ttreat = ttreat + 1;
  time++;
  
  LogicalVector host_reset = 
    (runif(pop_size, 0.0, 1.0) < 1-exp(-tau)) | (age > 1200);
  
  int num_reset = sum(host_reset);
  if(num_reset == 0) return;
  
  age[host_reset] = 0;
  WF[host_reset] = 0;
  WM[host_reset] = 0;
  Mf[host_reset] = 0;
  
  bite_risk[host_reset] = rgamma(num_reset,k,1/k);
  ttreat[host_reset] = 15*12;
  bednet[host_reset] = rbinom(num_reset, 1, bednet_coverage);
  double beta_shape1 = 
    treatment_coverage*(1-systematic_compliance)/systematic_compliance;
  double beta_shape2 = 
    (1-treatment_coverage)*(1-systematic_compliance)/systematic_compliance;
  treatment_prob[host_reset] = rbeta(num_reset, beta_shape1, beta_shape2);
  
  LogicalVector exclude = host_reset & (runif(pop_size, 0.0,1.0) < sys_ex);
  excluded[exclude] = 0;
}

void Data::importation() {
  LogicalVector importation_ind = 
    (runif(pop_size, 0.0, 1.0) < 1-exp(-1*importation_rate));
  if(is_false(any(importation_ind))) return;
  
  double mean_larvae = 10;
  int new_worms = floor(0.5 * lambda  * vth * psi_1 * psi_2 * s_2 *mean_larvae);
  //Rcout << "new worms = " << new_worms << "\n";
  WF[importation_ind] = new_worms;
  WM[importation_ind] = new_worms;
  Mf[importation_ind] = 0;
}

void Data::set_bednet(List bednet_data) {
  bool verbose = false;
  bednet_time = bednet_data["time"];
  bednet_coverage_vec = bednet_data["coverage"];
  if(verbose) {
    Rcout << "Bednet data set at times " << bednet_time
          << ", with coverage " << bednet_coverage_vec
          << "\n";
  }
}

void Data::set_MDA(List MDA_data) {
  bool verbose = false;
  MDA_time = MDA_data["time"];
  MDA_coverage = MDA_data["coverage"];
  if(verbose) {
    Rcout << "Setting MDA at times " << MDA_time
          << " with coverage " << MDA_coverage << "\n";
  }
}

void Data::set_parameter(std::string parameter, double value) {
  if(parameter == "k") {
    k = value;
    bite_risk = rgamma(pop_size, k, 1/k);
    return;
  }
  
  if(parameter == "vth") {
    vth = value;
    return;
  }
  
  if(parameter == "r") {
    r = value;
    return;
  }
  
  if(parameter == "alpha") {
    alpha = value;
    return;
  }
  
  if(parameter == "importation_rate") {
    importation_rate = value;
    return;
  }
  
  
  

  stop("Unknown parameter");
}

void Data::run_MDA() {
  if(!is_true(any(MDA_time == time))) return;
  bool verbose = false;
  if(verbose) {
    Rcout << "Running MDA at time t = " << time << "\n";
  }
  for(int i = 0; i < MDA_time.length(); i++) {
    if(time != MDA_time[i]) continue;
    if(treatment_coverage == MDA_coverage[i]) continue;
    
    // Redraw treatment probabilities
    treatment_coverage = MDA_coverage[i];
    double beta_shape1 = 
      treatment_coverage*(1-systematic_compliance)/systematic_compliance;
    double beta_shape2 = 
      (1-treatment_coverage)*(1-systematic_compliance)/systematic_compliance;
    treatment_prob = rbeta(pop_size, beta_shape1, beta_shape2);
  }
  NumericVector U = runif(pop_size, 0.0, 1.0);
  LogicalVector treat = (U < treatment_prob);

  for(int i = 0; i < pop_size; i++) {
    if(excluded.at(i) == 1) continue;
    if(treat.at(i) == 0) continue;
    if(age.at(i) < 5*12) continue;

    WM.at(i) = int(WM.at(i)*(1 - prop_worm_killed));
    WF.at(i) = int(WF.at(i)*(1 - prop_worm_killed));
    Mf.at(i) *= 1 - prop_mf_killed;

    ttreat.at(i) = 0;
  }
  larvae_dynamics();
}


void Data::update_bednet() {
  bool verbose = false;
  if(!is_true(any(bednet_time == time))) return;
  
  for(int i = 0; i < bednet_time.length(); i++) {
    int current_bednet_time = bednet_time.at(i);
    if(current_bednet_time != time) continue;
    
    bednet_coverage = bednet_coverage_vec.at(i);
    bednet = rbinom(pop_size, 1, bednet_coverage);
    if(verbose) {
      Rcout << "At time " << time 
            << ", updating bednets with new coverage = " << bednet_coverage
            << " and bednet = " << bednet << "\n";
    }
    return;
  }
}

void Data::burn_in() {
  for(int i = 0; i < 1200; i++) {
    timestep(1.0);
    time = 0;
  }
  int len = worm_prev.length();
  worm_prev = worm_prev[len-1];
  mf_prev = mf_prev[len-1];
  worm_prev67 = worm_prev67[len-1];
  mf_prev67 = mf_prev67[len-1];
  pop67 = pop67[len-1];
}

List Data::save() {
  List out;
  out["time"] = time;
  out["pop_size"] = pop_size;
  out["tau"] = tau;
  out["bite_risk"] = bite_risk;
  out["bednet"] = bednet;
  out["bednet_coverage"] = bednet_coverage;
  out["bednet_reduction"] = bednet_reduction;
  out["k"] = k;
  out["importation_rate"] = importation_rate;
  out["lambda"] = lambda;
  out["g"] = g;
  out["larvae"] = larvae;
  out["sigma"] = sigma;
  out["kappa"] = kappa;
  out["r"] = r;
  out["psi_1"]= psi_1;
  out["psi_2"] = psi_2;
  out["s_2"] = s_2;
  out["mu"] = mu;
  out["gamma"] = gamma;
  out["alpha"] = alpha;
  out["vth"] = vth;
  out["WF"] = WF;
  out["WM"] = WM;
  out["Mf"] = Mf;
  out["age"] = age;
  out["treatment_prob"] = treatment_prob;
  out["ttreat"] = ttreat;
  out["excluded"] = excluded;
  out["treatment_coverage"] = treatment_coverage;
  out["systematic_compliance"] = systematic_compliance;
  out["fecundity"] = fecundity;
  out["sys_ex"] = sys_ex;
  out["prop_mf_killed"] = prop_mf_killed;
  out["prop_worm_killed"] = prop_worm_killed;
  out["MDA_time"] = MDA_time;
  out["MDA_coverage"] = MDA_coverage;
  out["worm_prev"] = worm_prev;
  out["mf_prev"] = mf_prev;
  out["worm_prev67"] = worm_prev67;
  out["mf_prev67"] = mf_prev67;
  out["pop67"] = pop67;
  return clone(out);
}

void Data::load(const List input_data_raw) {
  List input_data = clone(input_data_raw);
  time = input_data["time"];
  pop_size = input_data["pop_size"];
  tau = input_data["tau"];
  bite_risk = input_data["bite_risk"];
  bednet= input_data["bednet"];
  bednet_coverage = input_data["bednet_coverage"];
  bednet_reduction = input_data["bednet_reduction"];
  k = input_data["k"];
  importation_rate = input_data["importation_rate"];
  lambda = input_data["lambda"];
  g = input_data["g"];
  larvae = input_data["larvae"];
  sigma = input_data["sigma"];
  kappa = input_data["kappa"];
  r = input_data["r"];
  psi_1 = input_data["psi_1"];
  psi_2 = input_data["psi_2"];
  s_2 = input_data["s_2"];
  mu = input_data["mu"];
  gamma = input_data["gamma"];
  alpha = input_data["alpha"];
  vth = input_data["vth"];
  WF = input_data["WF"];
  WM = input_data["WM"];
  Mf = input_data["Mf"];
  age = input_data["age"];
  treatment_prob = input_data["treatment_prob"];
  ttreat = input_data["ttreat"];
  excluded = input_data["excluded"];
  treatment_coverage = input_data["treatment_coverage"];
  systematic_compliance = input_data["systematic_compliance"];
  fecundity = input_data["fecundity"];
  sys_ex = input_data["sys_ex"];
  prop_mf_killed = input_data["prop_mf_killed"];
  prop_worm_killed = input_data["prop_worm_killed"];
  MDA_time = input_data["MDA_time"];
  MDA_coverage = input_data["MDA_coverage"];
  worm_prev = input_data["worm_prev"];
  mf_prev = input_data["mf_prev"];
  worm_prev67 = input_data["worm_prev67"];
  mf_prev67 = input_data["mf_prev67"];
  pop67 = input_data["pop67"];
}

