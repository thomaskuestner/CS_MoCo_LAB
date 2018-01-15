// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_PARSER_H
#define SP_PARSER_H
#include <map>

namespace sp
{
    ///
    /// @defgroup parser Parser
    /// \brief Parameter file parser functions.
    /// @{

    ///
    /// \brief A parser class.
    ///
    /// Implements parsing from text file for different types
    ///
    class parser
    {
        private:
            std::map<std::string, std::string> par_map; ///< Map structure to store parameter and value

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Check if key is valid.
            /// @returns TRUE if key is valid FALSE otherwise
            /// @param key Key string
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool valid_key(const std::string& key)
            {
                if(par_map.find(key) == par_map.end() )
                {
                    std::cout << "SigPack: Parameter "+ key + " not found!" <<std::endl;
                    return false;
                }
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Parse complex value.
            ///
            /// String may be in format e.g. "1+1i","1+1j" or entirely real or imaginary "1", "1i" or "1j"
            /// @returns Complex value
            /// @param str Complex notation string
            ////////////////////////////////////////////////////////////////////////////////////////////
            std::complex<double> parse_cx(std::string str)
            {
                double re,im;
                char i_ch;
                std::stringstream iss(str);

                // Parse full ...
                if(iss >> re >> im >> i_ch && (i_ch=='i'|| i_ch=='j')) return std::complex<double>(re,im);

                // ... or only imag
                iss.clear();
                iss.seekg(0,iss.beg);
                if(iss >> im >> i_ch && (i_ch=='i'|| i_ch=='j')) return std::complex<double>(0.0,im);

                // .. or only real
                iss.clear();
                iss.seekg(0,iss.beg);
                if(iss >> re) return std::complex<double>(re,0.0);

                // ... otherwise
                err_handler("Could not parse complex number!");
            }

        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ///
            /// Opens parameter file and puts the key and value in the map structure
            /// @param fname Parameter file name
            ////////////////////////////////////////////////////////////////////////////////////////////
            parser(const std::string& fname)
            {
                std::ifstream fh;
                std::string line;
                size_t mark = 0;

                // Clear parameter map
                par_map.clear();

                // Open file
                fh.open(fname.c_str());
                if (!fh)
                {
                    wrn_handler("Could not find " + fname);
                }
                else
                {
                    // Parse
                    while (std::getline(fh, line))
                    {
                        std::string keyS="";
                        std::string dataS="";

                        // Skip empty lines
                        if (line.empty())
                            continue;

                        // Skip lines with only whitespace
                        if(line.find_first_not_of("\t ")==std::string::npos)
                            continue;

                        // Remove comment
                        mark = line.find("%");
                        if(mark!=std::string::npos)
                            line.erase(mark,line.length());

                        // Do we have a '='
                        mark = line.find("=");
                        if(mark!=std::string::npos)
                        {
                            // Find key
                            keyS = line.substr(line.find_first_not_of("\t "),mark-line.find_first_not_of("\t "));
                            keyS = keyS.substr(0,keyS.find_last_not_of("\t ")+1);

                            // Find data
                            dataS = line.substr(mark+1,line.length());
                            dataS = dataS.substr(0,dataS.find_last_not_of("\t ")+1);

                            // Do we have a string
                            mark = dataS.find("\"");
                            if(mark!=std::string::npos)
                            {
                                dataS = dataS.substr(mark+1,dataS.length());
                                dataS = dataS.substr(0,dataS.find_last_of("\""));
                            }
                            // Do we have a vector/matrix
                            mark = dataS.find("[");
                            if(mark!=std::string::npos)
                            {
                                dataS = dataS.substr(mark+1,dataS.length());
                                dataS = dataS.substr(0,dataS.find_last_of("]"));
                            }

                            // Insert to map
                            par_map.insert(std::pair<std::string, std::string>(keyS, dataS));
                        }
                    }

                    // Close file
                    fh.close();
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~parser() {}

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Generic type get function.
            /// @returns Value of type T, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            T  getParam(const std::string key,const T def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default "<< def_val  <<std::endl;
                    return def_val;
                }
                std::istringstream iss(par_map.find(key)->second);
                T out;
                iss >> out;
                return out;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief String type get function.
            /// @returns String value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            std::string getString(const std::string key,const std::string def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default " << def_val <<std::endl;
                    return def_val;
                }
                return par_map.find(key)->second;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Col type get function.
            /// @returns Col value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            arma::Col<T> getCol(const std::string key,const arma::Col<T> def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string row,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),';')+1);
                arma::Col<T> x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, row, ';');
                    std::stringstream iss(row);
                    iss >> x(k);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief cx_vec type get function.
            /// @returns cx_vec value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_vec getCxCol(const std::string key,const arma::cx_vec def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string row,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),';')+1);
                arma::cx_vec x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, row, ';');
                    x(k) = parse_cx(row);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Row type get function.
            /// @returns Row value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            arma::Row<T> getRow(const std::string key,const arma::Row<T> def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string col,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
                arma::Row<T> x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, col, ',');
                    std::stringstream iss(col);
                    iss >> x(k);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief cx_rowvec type get function.
            /// @returns cx_rowvec value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_rowvec getCxRow(const std::string key,const arma::cx_rowvec def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string col,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
                arma::cx_rowvec x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, col, ',');
                    x(k) = parse_cx(col);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Mat type get function.
            /// @returns Mat value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            arma::Mat<T> getMat(const std::string key,const arma::Mat<T> def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }
                std::string full_str,row,col;
                std::istringstream iss_full;

                full_str = par_map.find(key)->second;
                int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

                iss_full.str(full_str);
                std::getline(iss_full, row, ';');
                int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

                arma::Mat<T> x(R,C);

                iss_full.seekg(0,iss_full.beg);
                for(int r=0; r<R; r++)
                {
                    std::getline(iss_full, row, ';');
                    std::istringstream iss_row(row);
                    for(int c=0; c<C; c++)
                    {
                        std::getline(iss_row, col, ',');
                        std::istringstream iss_col(col);
                        iss_col >> x(r,c);
                    }
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief cx_mat type get function.
            /// @returns cx_mat value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_mat getCxMat(const std::string key,const arma::cx_mat def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }
                std::string full_str,row,col;
                std::istringstream iss_full;

                full_str = par_map.find(key)->second;
                int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

                iss_full.str(full_str);
                std::getline(iss_full, row, ';');
                int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

                arma::cx_mat x(R,C);

                iss_full.seekg(0,iss_full.beg);
                for(int r=0; r<R; r++)
                {
                    std::getline(iss_full, row, ';');
                    std::istringstream iss_row(row);
                    for(int c=0; c<C; c++)
                    {
                        std::getline(iss_row, col, ',');
                        x(r,c)=parse_cx(col);
                    }
                }
                return x;
            }

    }; // end Class
    /// @}

} // end namespace
#endif
