// #include <RcppArmadillo.h>
// #include "Utils.h"
// 
// static double const log2pi = std::log(2.0 * M_PI);
// 
// class Parameters;
// class Data;
// class Transformations;
// 
// 
// class Parameters
// {
//   public:
//     arma::mat beta_map, beta_C, beta_chol;
//     arma::vec beta_g, beta_w, beta_mu;
//     arma::cube lambda_map, lambda_C, lambda_chol, lambda_old;
//     arma::mat lambda_g, lambda_w, lambda_mu, lambda_cov;
//     arma::cube eta_prec, eta_var;
//     arma::mat eta_map, eta_chol;
//     arma::vec eta_g, eta_w;
//     arma::vec prec_map;
//     double prec_a, prec_b;
//     arma::vec tau1;
//     arma::vec tau2;
//     Parameters(Data&);
//     void update_beta_map(Data&, Transformations&);
//     void update_lambda_map(Data&, Transformations&);
//     void update_eta_em(Data&, Transformations&);
//     void update_prec_map(Data&, Transformations&);
//     void update_tau1(Data&, Transformations&);
//     void update_tau2(Data&, Transformations&);
// };
// 
// class Transformations
// {
//   public:
//     arma::mat BtB, fit, lambda_tilde, BtY;
//     arma::cube ZetaZ;
//     arma::mat BlkDiagMeanPenalties;
//     arma::cube BlkDiagVarPenalties;
//     void BuildBlkDiagBeta(Data&, Parameters&);
//     void BuildBlkDiagLambda(Data&, Parameters&);
//     Transformations(Data&);
//     Transformations(Data&, Parameters&);
// };
// 
// class Data
// {
//   public:
//     arma::uword D1, D2, p, K, n, n_unique_mean, n_unique_var;
//     arma::mat Y, X, Z, B;
//     arma::uvec missing, missing_sub, missing_time;
//     arma::field<arma::mat> MeanPenalties;
//     arma::field<arma::mat> VarPenalties;
//     arma::uvec MeanIndices;
//     arma::uvec VarIndices;
//     Data(arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::uword,
//          arma::field<arma::mat>&,
//          arma::field<arma::mat>&,
//          arma::uvec&,
//          arma::uvec&);
//     void impute_response(Transformations&);
// };
// 
// Parameters::Parameters(Data& dat) {
//   beta_map = arma::mat(dat.p, dat.D1, arma::fill::randn);
//   beta_C = arma::mat(dat.p * dat.D1, dat.p * dat.D1);
//   beta_chol = arma::mat(dat.p * dat.D1, dat.p * dat.D1);
//   beta_g = arma::vec(dat.p * dat.D1);
//   beta_w = arma::vec(dat.p * dat.D1);
//   beta_mu = arma::vec(dat.p * dat.D1);
//   
//   lambda_map = arma::cube(dat.p, dat.D2, dat.K, arma::fill::randn);
//   lambda_C = arma::cube(dat.p * dat.D2, dat.p * dat.D2, dat.K);
//   lambda_chol = arma::cube(dat.p * dat.D2, dat.p * dat.D2, dat.K);
//   lambda_g = arma::mat(dat.p * dat.D2, dat.K);
//   lambda_w = arma::mat(dat.p * dat.D2, dat.K);
//   lambda_mu = arma::mat(dat.p * dat.D2, dat.K);
//   lambda_old = arma::cube(dat.p, dat.D2, dat.K);
//   eta_prec = arma::cube(dat.K, dat.K, dat.n);
//   eta_var = arma::cube(dat.K, dat.K, dat.n, arma::fill::ones);
//   eta_map = arma::mat(dat.n, dat.K, arma::fill::randn);
//   eta_chol = arma::mat(dat.K, dat.K);
//   eta_g = arma::vec(dat.K);
//   eta_w = arma::vec(dat.K);
//   
//   prec_map = arma::vec(dat.n, arma::fill::ones);
//   prec_map.fill(1);
//   prec_a = .0001;
//   prec_b = .0001;
//   tau1 = arma::vec(dat.MeanPenalties.n_elem, arma::fill::zeros);
//   tau2 = arma::vec(dat.VarPenalties.n_elem, arma::fill::zeros);
// }
// 
// Transformations::Transformations(Data& dat, Parameters& pars) {
//   BtY = dat.B.t() * dat.Y.t();
//   BtB = dat.B.t() * dat.B;
//   // initialize fit
//   fit = dat.X * pars.beta_map.t() * dat.B.t();
//   for (arma::uword k = 0; k < dat.K; k++) {
//     fit = fit + arma::diagmat(pars.eta_map.col(k)) * 
//       dat.Z * pars.lambda_map.slice(k).t() * dat.B.t();
//   }
//   // initialize lambda_tilde
//   lambda_tilde = arma::mat(dat.p, dat.K);
//   // initialize ZetaZ
//   arma::vec eta_tube;
//   ZetaZ = arma::cube(dat.D2, dat.D2, dat.K);
//   for (arma::uword k = 0; k < dat.K; k++) {
//     eta_tube = pars.eta_var.tube(k, k);
//     ZetaZ.slice(k) = dat.Z.t() * arma::diagmat(pars.prec_map) *
//       arma::diagmat(arma::square(pars.eta_map.col(k)) +
//       eta_tube) * dat.Z;
//   }
//   
//   BlkDiagMeanPenalties = arma::mat(
//     dat.p * dat.D1, dat.p * dat.D1, arma::fill::zeros);
//   BlkDiagVarPenalties = arma::cube(
//     dat.p * dat.D2, dat.p * dat.D2, dat.K, arma::fill::zeros);
// }
// Data::Data(arma::mat& y, arma::mat& x, arma::mat& z, arma::mat& b,
//            arma::uword dim, arma::field<arma::mat>& mp,
//            arma::field<arma::mat>& vp,
//            arma::uvec& mi,
//            arma::uvec& vi) {
//   D1 = x.n_cols;
//   D2 = z.n_cols;
//   p = b.n_cols;
//   K = dim;
//   n = y.n_rows;
//   Y = y;
//   X = x;
//   Z = z;
//   B = b;
//   Y = y;
//   missing = arma::find_nonfinite(Y);
//   missing_sub = armadillo_modulus3(missing, Y.n_rows);
//   missing_time = arma::floor(missing / Y.n_rows);
//   Y.elem(missing).fill(0);
//   MeanPenalties = mp;
//   VarPenalties = vp;
//   MeanIndices = mi;
//   VarIndices = vi;
//   
// }
// 
// void Parameters::update_beta_map(Data& dat, Transformations& transf) {
//   beta_C = arma::kron(dat.X.t() * arma::diagmat(prec_map) * dat.X,
//                       transf.BtB) + transf.BlkDiagMeanPenalties;
//   beta_chol = arma::chol(beta_C, "lower");
//   beta_g = arma::vectorise(dat.B.t() * (dat.Y.t() - transf.fit.t() + 
//     dat.B * beta_map * dat.X.t()) * arma::diagmat(prec_map) * dat.X);
//   beta_w = arma::solve(arma::trimatl(beta_chol), beta_g);
//   beta_mu = arma::solve(arma::trimatu(beta_chol.t()), beta_w);
//   // update fit
//   transf.fit = transf.fit - dat.X * beta_map.t() * dat.B.t();
//   beta_map = arma::reshape(beta_mu, dat.p, dat.D1);
//   transf.fit = transf.fit + dat.X * beta_map.t() * dat.B.t();
// }
// 
// //void Parameters::update_lambda_map_whole(Data& dat, Transformations& transf) {
//   
// //}
// void Parameters::update_lambda_map(Data& dat, Transformations& transf) {
//   lambda_old = lambda_map;
//   arma::vec eta_cov;
//   for (arma::uword k = 0; k < dat.K; k++) {
//     lambda_C.slice(k) = arma::kron(transf.ZetaZ.slice(k), transf.BtB) +
//       transf.BlkDiagVarPenalties.slice(k);
//     if (!lambda_C.slice(k).is_symmetric()) {
//       lambda_C.slice(k) = arma::symmatu(lambda_C.slice(k));
//     }
//     lambda_chol.slice(k) = arma::chol(lambda_C.slice(k), "lower");
//     lambda_g.col(k) = arma::vectorise(dat.B.t() * (dat.Y.t() - dat.B *
//       beta_map * dat.X.t()) * arma::diagmat(eta_map.col(k)) * 
//       arma::diagmat(prec_map) * dat.Z);
//     for (arma::uword kp = 0; kp < dat.K; kp++) {
//       if (k != kp) {
//         eta_cov = eta_var.tube(k, kp);
//         
//         lambda_g.col(k) = lambda_g.col(k) - 
//           arma::vectorise(transf.BtB * lambda_map.slice(kp) * 
//           dat.Z.t() * (arma::diagmat(eta_map.col(kp)) *
//           arma::diagmat(eta_map.col(k)) +
//           arma::diagmat(eta_cov)) * arma::diagmat(prec_map) * dat.Z);
//       }
//     }/*
//     lambda_g.col(k) = arma::vectorise(dat.B.t() * (((dat.Y.t() -
//       transf.fit.t() + dat.B * lambda_map.slice(k) * 
//       dat.Z.t() * arma::diagmat(eta_map.col(k))) * arma::diagmat(prec_map)) * 
//       arma::diagmat(eta_map.col(k)) - 1.0 / 2.0 * (lambda_cov)) * dat.Z);*/
//     lambda_w.col(k) = arma::solve(arma::trimatl(lambda_chol.slice(k)), lambda_g.col(k));
//     lambda_mu.col(k) = arma::solve(arma::trimatu(lambda_chol.slice(k).t()), lambda_w.col(k));
//     // update fit
//     transf.fit = transf.fit - arma::diagmat(eta_map.col(k)) * dat.Z * 
//       lambda_map.slice(k).t() * dat.B.t();
//     lambda_map.slice(k) = arma::reshape(lambda_mu.col(k), dat.p, dat.D2);
//     transf.fit = transf.fit + arma::diagmat(eta_map.col(k)) * dat.Z * 
//       lambda_map.slice(k).t() * dat.B.t();
//     
//   }
// }
// 
// void Parameters::update_eta_em(Data& dat, Transformations& transf) {
//   for (arma::uword i = 0; i < dat.n; i++) {
//     for (arma::uword k = 0; k < dat.K; k++) {
//       transf.lambda_tilde.col(k) = lambda_map.slice(k) * dat.Z.row(i).t();
//     }
//     eta_prec.slice(i) = prec_map(i) * transf.lambda_tilde.t() * transf.BtB
//       * transf.lambda_tilde + arma::eye<arma::mat>(dat.K, dat.K);
//     if (!eta_prec.slice(i).is_symmetric()) {
//       eta_prec.slice(i) = arma::symmatu(eta_prec.slice(i));
//     }
//     eta_chol = arma::chol(eta_prec.slice(i), "lower");
//     eta_g = transf.lambda_tilde.t() * dat.B.t() *  
//       prec_map(i) * (dat.Y.row(i).t() - dat.B * beta_map * dat.X.row(i).t());
//     eta_w = arma::solve(arma::trimatl(eta_chol), eta_g);
//     eta_map.row(i) = arma::solve(arma::trimatu(eta_chol.t()), eta_w).t();
//     eta_var.slice(i) = arma::inv_sympd(eta_prec.slice(i));
//   }
//   
//   
//   // updated fit
//   transf.fit = dat.X * beta_map.t() * dat.B.t();
//   for (arma::uword k = 0; k < dat.K; k++) {
//     transf.fit = transf.fit + arma::diagmat(eta_map.col(k)) * 
//       dat.Z * lambda_map.slice(k).t() * dat.B.t();
//   }
//   
//   // updated ZetaZ
//   arma::vec eta_tube;
//   for (arma::uword k = 0; k < dat.K; k++) {
//     eta_tube = eta_var.tube(k, k);
//     transf.ZetaZ.slice(k) = dat.Z.t() * arma::diagmat(prec_map) *
//       arma::diagmat(arma::square(eta_map.col(k)) +
//       eta_tube) * dat.Z;
//   }
// }
// void Parameters::update_prec_map(Data& dat, Transformations& transf) {
//   double numerator = prec_a + double(dat.B.n_rows) / 2 - 1;
//   //double numerator = double(dat.B.n_rows) / 2;
//   double denominator;
//   arma::rowvec z;
//   for (arma::uword i = 0; i < dat.n; i++) {
//     for (arma::uword k = 0; k < dat.K; k++) {
//       transf.lambda_tilde.col(k) = lambda_map.slice(k) * dat.Z.row(i).t();
//     }
//     z = dat.Y.row(i) - transf.fit.row(i);
//     //denominator = arma::dot(z, z);
//     //denominator = denominator + 
//     //denominator = denominator + arma::as_scalar(eta_map.row(i) *
//     //  transf.lambda_tilde.t() * transf.BtB * transf.lambda_tilde * eta_map.row(i).t()) + 1.0 / 2.0 *
//     //  arma::trace(transf.lambda_tilde.t() * transf.BtB * transf.lambda_tilde * eta_var.slice(i));
//     z = dat.Y.row(i) - arma::trans(dat.B * beta_map * dat.X.row(i).t());
//     denominator = arma::dot(z,z) - 2 * arma::as_scalar(z * dat.B * transf.lambda_tilde * eta_map.row(i).t());
//     denominator = denominator + arma::as_scalar(eta_map.row(i) *
//       transf.lambda_tilde.t() * transf.BtB * transf.lambda_tilde * eta_map.row(i).t()) +  
//       arma::trace(transf.lambda_tilde.t() * transf.BtB * transf.lambda_tilde * eta_var.slice(i));
//     denominator = prec_b + 0.5 * denominator;
//     //denominator = 0.5 * denominator;
//     prec_map(i) = numerator / denominator;
//   }
// }
// 
// void Transformations::BuildBlkDiagBeta(Data& dat, Parameters& pars) {
//   arma::uword old_index = 0;
//   arma::uword start;
//   arma::uword end = -1;
//   for (arma::uword d = 0; d < dat.MeanIndices.n_elem; d++) {
//     if (dat.MeanIndices(d) == old_index) {
//       BlkDiagMeanPenalties.submat(
//         start, start, end, end) =
//           pars.tau1(d) * dat.MeanPenalties(d) + BlkDiagMeanPenalties.submat(
//             start, start, end, end);
//     } else {
//       start = end + 1;
//       end = start + dat.MeanPenalties(d).n_rows - 1;
//       BlkDiagMeanPenalties.submat(
//         start, start, end, end) =
//           pars.tau1(d) * dat.MeanPenalties(d);
//       
//     }
//     
//     old_index = dat.MeanIndices(d);
//   }
//   
// }
// void Transformations::BuildBlkDiagLambda(Data& dat, Parameters& pars) {
//   for (arma::uword k = 0; k < dat.K; k++) {
//     arma::uword old_index = 0;
//     arma::uword start;
//     arma::uword end = -1;
//     for (arma::uword d = 0; d < dat.VarIndices.n_elem; d++) {
//       if (dat.VarIndices(d) == old_index) {
//         BlkDiagVarPenalties.slice(k).submat(
//           start, start, end, end) =
//             pars.tau2(d) * dat.VarPenalties(d) + BlkDiagVarPenalties.slice(k).submat(
//                 start, start, end, end);
//       } else {
//         start = end + 1;
//         end = start + dat.VarPenalties(d).n_rows - 1;
//         BlkDiagVarPenalties.slice(k).submat(
//           start, start, end, end) =
//             pars.tau2(d) * dat.VarPenalties(d);
//         
//       }
//       
//       old_index = dat.VarIndices(d);
//     }
//   }
//   
// }
// void Parameters::update_tau1(Data& dat, Transformations& transf) {
//   
//   
//   double a = 1;
//   double b = .005;
//   double update_a = 0, update_b = 0;
//   arma::uword start = 0;
//   arma::uword end = static_cast<double>(dat.MeanPenalties(0).n_rows) /
//     static_cast<double>(dat.p) - 1;
//   arma::uword num_field_elements = dat.MeanPenalties.n_elem;
//   arma::uword old_index = 1;
//   
//   for(arma::uword i = 0; i < num_field_elements; i++){
//     
//     if(dat.MeanIndices(i) != old_index){
//       start = end + 1;
//       end = end + dat.MeanPenalties(i).n_rows / dat.p;
//     } 
//     
//     update_a = dat.MeanPenalties(i).n_rows;
//     update_b = 1.0 / 2.0 * 
//       arma::as_scalar(arma::vectorise(beta_map.cols(start, end)).t() *
//       dat.MeanPenalties(i) *
//       arma::vectorise(beta_map.cols(start, end)));
//     tau1(i) = (a + update_a / 2.0 - 1) / (b + update_b);
//     old_index = dat.MeanIndices(i);
//   }
//   transf.BuildBlkDiagBeta(dat, *this);
//   
// }
// 
// void Parameters::update_tau2(Data& dat, Transformations& transf) {
//   double a = 1;
//   double b = .005;
//   double update_a = 0, update_b = 0;
//   arma::uword start = 0;
//   arma::uword end = static_cast<double>(dat.VarPenalties(0).n_rows) /
//     static_cast<double>(dat.p) - 1;
//   arma::uword num_field_elements = dat.VarPenalties.n_elem;
//   arma::uword old_index = 1;
//   
//   for(arma::uword i = 0; i < num_field_elements; i++){
//     
//     if (dat.VarIndices(i) != old_index) {
//       start = end + 1;
//       end = end + dat.VarPenalties(i).n_rows / dat.p;
//     } 
// 
//     for (arma::uword k = 0; k < dat.K; k++) {
//       update_a = update_a + dat.VarPenalties(i).n_rows;
//       update_b = update_b + 
//         arma::as_scalar(arma::vectorise(lambda_map.slice(k).cols(start, end)).t() *
//         dat.VarPenalties(i) * arma::vectorise(lambda_map.slice(k).cols(start, end)));
//     }
//     tau2(i) = (a + update_a / 2.0 - 1) / (b + 0.5 * update_b);
//     update_a = 0;
//     update_b = 0;
//     old_index = dat.VarIndices(i);
//     
//   }
//   transf.BuildBlkDiagLambda(dat, *this);
// }
// 
// void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
//   arma::uword const n = trimat.n_cols;
//   
//   for(unsigned j = n; j-- > 0;){
//     double tmp(0.);
//     for(unsigned i = 0; i <= j; ++i)
//       tmp += trimat.at(i, j) * x[i];
//     x[j] = tmp;
//   }
// }
// 
// // [[Rcpp::export]]
// arma::vec dmvnrm_arma_fast(arma::mat const &x,  
//                            arma::rowvec const &mean,  
//                            arma::mat const &sigma, 
//                            bool const logd = false) { 
//   using arma::uword;
//   uword const n = x.n_rows, 
//     xdim = x.n_cols;
//   arma::vec out(n);
//   arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
//   double const rootisum = arma::sum(log(rooti.diag())), 
//     constants = -(double)xdim/2.0 * log2pi, 
//     other_terms = rootisum + constants;
//   
//   arma::rowvec z;
//   for (uword i = 0; i < n; i++) {
//     z = (x.row(i) - mean);
//     inplace_tri_mat_mult(z, rooti);
//     out(i) = other_terms - 0.5 * arma::dot(z, z);     
//   }  
//   
//   if (logd)
//     return out;
//   return exp(out);
// }
// 
// 
// double compute_posterior_density(Parameters& pars, Data& dat,
//                                  Transformations& transf) {
//   double log_post = 0;
//   arma::rowvec mu;
//   arma::mat cov;
//   arma::mat identity_mat = arma::eye<arma::mat>(dat.B.n_rows, dat.B.n_rows);
//   arma::mat z;
//   arma::mat small_cov(dat.B.n_cols, dat.B.n_cols);
//   for (arma::uword i = 0; i < dat.n; i++) {
//     small_cov.zeros();
//     mu = arma::trans(dat.B * pars.beta_map * dat.X.row(i).t());
//     for (arma::uword k = 0; k < dat.K; k++) {
//       z = pars.lambda_map.slice(k) * dat.Z.row(i).t();
//       small_cov = small_cov + z * z.t();
//     }
//     cov = dat.B * small_cov * dat.B.t() + 1.0 /
//       pars.prec_map(i) * arma::diagmat(identity_mat);
//     log_post = log_post + arma::as_scalar(
//       dmvnrm_arma_fast(dat.Y.row(i), mu, cov, true));
//   }
//   for (arma::uword i = 0; i < dat.n; i++) {
//     log_post = log_post + R::dgamma(pars.prec_map(i), pars.prec_a, 1 / pars.prec_b, 1);
//   }
//   
//   arma::vec eigval;
//   arma::mat eigmat;
//   double log_det_beta = 0;
//   arma::eig_sym(eigval, eigmat, transf.BlkDiagMeanPenalties);
//   for (arma::uword i = 0; i < eigval.n_elem; i++) {
//     if (eigval(i) > 1e-6) {
//       log_det_beta = std::log(eigval(i)) + log_det_beta;
//     }
//   }
//   log_det_beta = .5 * (log_det_beta -
//     static_cast<double>(eigval.n_elem) * log2pi);
//   log_post = log_post + log_det_beta -.5 *
//     arma::as_scalar(arma::vectorise(pars.beta_map).t() *
//     transf.BlkDiagMeanPenalties * arma::vectorise(pars.beta_map));
//   
//   
//   for (arma::uword d = 0; d < dat.MeanIndices.n_elem; d++) {
//     log_post = log_post + R::dgamma(pars.tau1(d), 1, 1/.005, 1);
//   }
// 
//   
//   
// 
//   
//   arma::vec eigval_lambda;
//   arma::mat eigmat_lambda;
//   double log_det_lambda = 0;
//   arma::eig_sym(eigval_lambda, eigmat_lambda, transf.BlkDiagVarPenalties.slice(0));
//   for (arma::uword i = 0; i < eigval_lambda.n_elem; i++) {
//     if (eigval(i) > 1e-6) {
//       log_det_lambda = std::log(eigval_lambda(i)) + log_det_lambda;
//     }
//   }
//   log_det_lambda = .5 * (log_det_lambda -
//     static_cast<double>(eigval_lambda.n_elem) * log2pi);
//   for (arma::uword k = 0; k < dat.K; k++) {
//     log_post = log_post + log_det_lambda + (-1.0 / 2.0 * 
//       arma::as_scalar(arma::vectorise(pars.lambda_map.slice(k)).t() * 
//       transf.BlkDiagVarPenalties.slice(0) * arma::vectorise(pars.lambda_map.slice(k))));
//   }
//   
//   for (arma::uword d = 0; d < dat.VarIndices.n_elem; d++) {
//     log_post = log_post + R::dgamma(pars.tau2(d), 1, 1 / .005, 1);
//   }
//   
// 
//   return(log_post);
// }
// 
// // [[Rcpp::export]]
// Rcpp::List my_main(arma::mat& Y, arma::mat& X, arma::mat& Z,
//                    arma::mat& B, arma::uword K, arma::uword iter,
//                    arma::field<arma::mat>& MeanPenalties,
//                    arma::field<arma::mat>& VarPenalties,
//                    arma::uvec& MeanIndices,
//                    arma::uvec& VarIndices) {
//   double conv;
//   Data dat(Y, X, Z, B, K, MeanPenalties, VarPenalties, MeanIndices, VarIndices);
//   Rcpp::Rcout << "Initialized data" << std::endl; 
//   Parameters pars(dat);
//   Rcpp::Rcout << "Initialized parameters" << std::endl; 
//   Transformations transf(dat, pars);
//   Rcpp::Rcout << "Initialized transformations" << std::endl; 
//   arma::vec log_post(iter);
//   
//   transf.BuildBlkDiagBeta(dat, pars);
//   transf.BuildBlkDiagLambda(dat, pars);
//   for (arma::uword i = 0; i < iter; i++) {
//     Rcpp::Rcout << i << std::endl;
//     pars.update_eta_em(dat, transf);
//     pars.update_beta_map(dat, transf);
//     pars.update_lambda_map(dat, transf);
//     pars.update_prec_map(dat, transf);
//     pars.update_tau1(dat, transf);
//     pars.update_tau2(dat, transf);
//     //log_post(i) = compute_posterior_density(pars, dat, transf);
//     conv = arma::accu(arma::square(pars.lambda_map - pars.lambda_old)) /
//       arma::accu(arma::square((pars.lambda_old)));
//     Rcpp::Rcout << conv << std::endl;
//     //Rcpp::Rcout << log_post(i) <<  std::endl;
//   }
//   //log_post(0) = compute_posterior_density(pars, dat, transf);
//   return(Rcpp::List::create(Rcpp::Named("Beta", pars.beta_map),
//                             Rcpp::Named("Lambda", pars.lambda_map),
//                             Rcpp::Named("Precision", pars.prec_map),
//                             Rcpp::Named("Eta", pars.eta_map),
//                             Rcpp::Named("log_post", log_post),
//                             Rcpp::Named("fit", transf.fit),
//                             Rcpp::Named("Eta_var", pars.eta_var),
//                             Rcpp::Named("Eta_g", pars.eta_g),
//                             Rcpp::Named("tau1", pars.tau1),
//                             Rcpp::Named("tau2", pars.tau2),
//                             Rcpp::Named("BlkdDiagMean", transf.BlkDiagMeanPenalties),
//                             Rcpp::Named("BlkDiagVar", transf.BlkDiagVarPenalties)));
// }
// 
