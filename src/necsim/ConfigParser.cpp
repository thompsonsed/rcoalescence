//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
// 
/**
 * @author Sam Thompson
 * @date 31/08/2016
 * @file ConfigFileParser.cpp
 * @brief Contains implementation of the ConfigFileParser.h functions.
 * 
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 *
 */



#include "cpp17_includes.h"
#include "custom_exceptions.h"
#include "Logging.h"
#include "ConfigParser.h"

namespace necsim
{
    void importArgs(const unsigned int &argc, char* argv[], std::vector<string> &comargs)
    {
        for(unsigned int i = 0; i < argc; i++)
        {
            comargs.emplace_back(argv[i]);
        }
        // check size is correct
        if(comargs.size() != argc)
        {
            writeError("ERROR_MAIN_010: Incorrect command line parsing.");
        }
    }

    std::string SectionOption::getOption(std::string refval)
    {
        for(unsigned int i = 0; i < refs.size(); i++)
        {
            if(refs[i] == refval)
            {
                return (val[i]);
            }
        }
#ifdef DEBUG
        std::stringstream ss;
        ss << "Reference " << refval << " not found in keyoption." << std::endl;
        writeInfo(ss.str());
#endif
        return ("null");
    }

    std::ostream &operator<<(std::ostream &os, const SectionOption &k)
    {
        os << k.section << "\n" << k.val.size() << "\n" << k.refs.size() << "\n";
        for(const auto &i : k.val)
        {
            os << i << "\n";
        }
        for(const auto &ref : k.refs)
        {
            os << ref << "\n";
        }
        return os;
    }

    std::istream &operator>>(std::istream &is, SectionOption &k)
    {
        // os << m.num_rows<<" , "<<m.num_cols<<" , "<<std::endl;
        unsigned int valsize, refsize;
        is >> k.section >> valsize >> refsize;
        is.ignore();
        std::string tmp;
        for(unsigned int i = 0; i < valsize; i++)
        {
            getline(is, tmp);
            k.val.push_back(tmp);
        }
        for(unsigned int i = 0; i < refsize; i++)
        {
            getline(is, tmp);
            k.refs.push_back(tmp);
        }
        return is;
    }

    void ConfigParser::setConfig(const std::string &file, bool main, bool full_parse)
    {
        if(!configSet)
        {
            isMain = main;
            config_file = file;
            configSet = true;
            isFullParser = full_parse;
        }
        else
        {
            throw ConfigException("Attempt to set config file twice.");
        }
    }

    void ConfigParser::parseConfig()
    {
        std::ifstream is_file;
        if(!fs::exists(config_file))
        {
            std::stringstream ss;
            ss << "No config file found at " << config_file << ". Check file exists." << std::endl;
            throw ConfigException(ss.str());
        }
        try
        {
            is_file.open(config_file);
        }
        catch(...)
        {
            throw ConfigException("ERROR_CONF_004c: Could not open the config file. Check file exists and is readable.");
        }
        parseConfig(is_file);
        if(is_file.eof())
        {
            is_file.close();
        }
        else
        {
            throw ConfigException("ERROR_CONF_002: End of file not reached. Check input file formatting.");
        }
    }

    void ConfigParser::parseConfig(std::istream &istream)
    {
        if(!istream.fail() || !istream.good())
        {
            std::string line;
            // Get the first line of the file.
            while(getline(istream, line))
            {
                std::istringstream is_line(line);
                std::string key;
                std::string val;
                // Skip all whitespace
                is_line >> std::skipws;
                // start a new section
                if(line[0] == '[')
                {
                    SectionOption tempSections;
                    // get the section name
                    std::string section;
                    if(getline(is_line, section, ']'))
                    {
                        section = section.erase(0, 1);
                        tempSections.section = section;
                    }
                    // read each line
                    while(getline(istream, line))
                    {
                        // end the section when a new one starts.
                        if(line[0] == '[' || line.size() == 0)
                        {
                            break;
                        }
                        std::istringstream is_line2(line); // update the input-line stream
                        if(getline(is_line2, key, '='))
                        {

                            key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                            is_line2 >> std::skipws;
                        }
                        if(!is_line2)
                        {
                            throw ConfigException("ERROR_CONF_001: Read error in config file.");
                        }
                        if(getline(is_line2, val))
                        {
                            //						 	This line has been removed to allow for white spaces in file names and paths
                            //							val.erase(std::remove(val.begin(), val.end(), ' '), val.end());
                            while(val[0] == ' ')
                            {
                                val.erase(val.begin(), val.begin() + 1);
                            }
                        }
                        if(!is_line2)
                        {
                            throw ConfigException("ERROR_CONF_001: Read error in config file.");
                        }
                        tempSections.refs.push_back(key);
                        tempSections.val.push_back(val);
                    }
                    configs.push_back(tempSections);
                }
            }
        }
        else
        {
            throw ConfigException("ERROR_CONF_004b: Could not open the config file " + config_file
                                  + ". Check file exists and is readable.");
        }
    }

    std::vector<SectionOption> ConfigParser::getSectionOptions()
    {
        return configs;
    }

    void ConfigParser::setSectionOption(std::string section, std::string reference, std::string value)
    {
        SectionOption* section_option = nullptr;
        for(auto &option : configs)
        {
            if(option.section == section)
            {
                section_option = &option;
                break;
            }
        }
        if(!section_option)
        {
            SectionOption tmp;
            tmp.section = section;
            configs.emplace_back(tmp);
            section_option = &configs.back();
        }
        section_option->refs.emplace_back(reference);
        section_option->val.emplace_back(value);
    }

    SectionOption ConfigParser::operator[](unsigned long index)
    {
        return configs[index];
    }

    unsigned long ConfigParser::getSectionOptionsSize()
    {
        return configs.size();
    }

    std::vector<string> ConfigParser::getSections()
    {
        std::vector<string> toret;
        for(auto &config : configs)
        {
            toret.push_back(config.section);
        }
        return toret;
    }

    bool ConfigParser::hasSection(const std::string &sec)
    {
        for(auto &config : configs)
        {
            if(config.section == sec)
            {
                return true;
            }
        }
        return false;
    }

    std::vector<string> ConfigParser::getSectionValues(std::string sec)
    {
        for(auto &config : configs)
        {
            if(config.section == sec)
            {
                return config.val;
            }
        }
        throw ConfigException("Section not found in config file: " + sec);
    }

    std::string ConfigParser::getSectionOptions(std::string section, std::string ref)
    {
        for(auto &config : configs)
        {
            if(config.section == section)
            {
                for(unsigned int j = 0; j < config.refs.size(); j++)
                {
                    if(config.refs[j] == ref)
                    {
                        return config.val[j];
                    }
                }
            }
        }
#ifdef DEBUG
        writeWarning("No reference found for " + section + ", ");
#endif
        return "null";
    }

    std::string ConfigParser::getSectionOptions(std::string section, std::string ref, std::string def)
    {
        for(auto &config : configs)
        {
            if(config.section == section)
            {
                for(unsigned int j = 0; j < config.refs.size(); j++)
                {
                    if(config.refs[j] == ref)
                    {
                        return config.val[j];
                    }
                }
            }
        }
        return def;
    }

    int ConfigParser::importConfig(std::vector<string> &comargs)
    {
        // Check that the previous arguments have already been imported.
        if(isMain)
        {
            if(comargs.size() != 3)
            {
                throw ConfigException("ERROR_CONF_003: Number of command line arguments not correct before import.");
            }
        }
        std::ifstream is_file;
        try
        {
            is_file.open(config_file);
        }
        catch(...)
        {
            throw ConfigException("ERROR_CONF_004a: Could not open the config file. Check file exists and is readable.");
        }
        if(!is_file.fail())
        {
            std::string line;
            while(getline(is_file, line))
            {
                std::istringstream is_line(line);
                std::string key;
                is_line >> std::skipws;
                if(line[0] == '[')
                {
                    continue;
                }
                if(getline(is_line, key, '='))
                {
                    // Could implement proper data parsing based on the key object.
                    is_line >> std::skipws;
                    std::string value;
                    if(getline(is_line, value))
                    {
                        value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                        if(!is_line)
                        {
                            std::stringstream os;
                            os << value << std::endl;
                            writeWarning(os.str());
                            throw ConfigException("ERROR_CONF_001: Read error in config file.");
                        }
                        auto* tmp = new char[value.length() + 1];
                        strcpy(tmp, value.c_str());
                        comargs.emplace_back(tmp);
                    }
                }
            }
        }
        else
        {
            throw ConfigException("ERROR_CONF_004d: Could not open the config file. Check file exists and is readable.");
        }
        if(is_file.eof())
        {
            is_file.close();
        }
        else
        {
            throw ConfigException("ERROR_CONF_002: End of file not reached. Check input file formating.");
        }
        if(isMain)
        {
            // remove the file name from the command line arguments to maintain the std::vector format.
            comargs.erase(comargs.begin() + 2);
        }
        return static_cast<int>(comargs.size());
    }

    std::ostream &operator<<(std::ostream &os, const ConfigParser &c)
    {
        os << c.config_file << "\n" << c.configSet << "\n" << c.isMain << "\n" << c.isFullParser << "\n"
           << c.configs.size() << "\n";
        for(const auto &config : c.configs)
        {
            os << config;
        }
        return os;
    }

    std::istream &operator>>(std::istream &is, ConfigParser &c)
    {
        unsigned int configsize;
        is.ignore();
        getline(is, c.config_file);
        is >> c.configSet >> c.isMain >> c.isFullParser >> configsize;
        SectionOption tmpoption;
        if(configsize > 10000)
        {
            throw std::runtime_error("Config size extremely large, check file: " + std::to_string(configsize));
        }
        if(configsize > 0)
        {
            for(unsigned int i = 0; i < configsize; i++)
            {
                is >> tmpoption;
                c.configs.push_back(tmpoption);
            }
        }
        //		os << "end config" << std::endl;
        return is;
    }
}