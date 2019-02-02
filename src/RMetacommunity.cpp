////
//// Created by Sam Thompson on 28/01/2019.
////
//
//#include "RMetacommunity.h"
//
//RMetacommunity::RMetacommunity() : spec_sim_parameters(make_shared<SpecSimParameters>())
//{
//
//}
//
//void RMetacommunity::setup(const string &database, const bool &record_spatial, const string &sample_file,
//                           const string &use_fragments_in)
//{
//    spec_sim_parameters->setup(database, record_spatial, sample_file, use_fragments_in);
//}
//
//void RMetacommunity::addSpeciationRate(const double &spec_rate)
//{
//    spec_sim_parameters->addSpeciationRate(spec_rate);
//}
//
//void RMetacommunity::addTime(const double &time)
//{
//    spec_sim_parameters->addTime(time);
//}
//
//void RMetacommunity::addMetacommunityParameters(const unsigned long &metacommunity_size_in,
//                                                const double &metacommunity_speciation_rate_in,
//                                                const string &metacommunity_option_in,
//                                                const unsigned long &metacommunity_reference_in)
//{
//    spec_sim_parameters->addMetacommunityParameters(metacommunity_size_in, metacommunity_speciation_rate_in,
//                                                    metacommunity_option_in, metacommunity_reference_in);
//}
//
//void RMetacommunity::addProtractedParameters(const double &proc_spec_min, const double &proc_spec_max)
//{
//    spec_sim_parameters->addProtractedParameters(proc_spec_min, proc_spec_max);
//}
//
//void RMetacommunity::apply()
//{
//    applyNoOutput(spec_sim_parameters);
//    output();
//}
