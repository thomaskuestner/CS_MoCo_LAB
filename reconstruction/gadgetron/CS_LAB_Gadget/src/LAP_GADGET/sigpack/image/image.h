#ifndef IMAGE_H_INCLUDED
#define IMAGE_H_INCLUDED

#include <string>
#include <fstream>
#include <bitset>

using namespace std;
namespace sp
{
    ///
    /// @defgroup image Image
    /// \brief Image functions.
    /// @{

    ///
    /// \brief Portable anymap format class.
    ///
    /// Implements portable anymap image functions
    /// Supports .pbm, .pgm and .ppm plain and raw
    ///

    class PNM
    {
        private:
            std::ifstream ifs;          ///< Input stream handle
            std::ofstream ofs;          ///< Output stream handle
            arma::uword cols;                   ///< Nr of columns in image
            arma::uword rows;                   ///< Nr of rows in image
            int maxval;                 ///< Maximum pixel value in image
        public:
            enum imtype { NOTUSED, PBM_A, PGM_A, PPM_A, PBM_B, PGM_B, PPM_B } type; ///< Image format

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            PNM()
            {
                type   = NOTUSED;
                cols   = 0;
                rows   = 0;
                maxval = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~PNM() {}

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Clears the internal variables.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void clear(void)
            {
                type   = NOTUSED;
                cols   = 0;
                rows   = 0;
                maxval = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Reads the .pnm header.
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////
            void read_header()
            {
                string str;
                while (ifs >> str ) // Read until we have maxval
                {
                    // Remove comments
                    size_t mark = str.find_first_of("#");
                    if ( mark!= string::npos)
                    {
                        ifs.ignore(256, '\n');
                        str.erase(mark, 1);

                        if (str.empty()) continue;
                    }

                    if (type == NOTUSED)
                    {
                        type = static_cast<imtype>(str.at(1)-48);  // Conv char to int
                        if(str.at(0)!='P' || type>PPM_B) err_handler("Wrong type!");
                    }
                    else if (cols   == 0)
                    {
                        cols = atoi(str.c_str());
                    }
                    else if (rows   == 0)
                    {
                        rows = atoi(str.c_str());
                        if(type==PBM_A || type==PBM_B)
                        {
                            maxval = 1;  // Set maxval=1 for .PBM types
                            break;
                        }
                    }
                    else if (maxval == 0)
                    {
                        maxval = atoi(str.c_str());
                        break;
                    }
                }
                ifs.ignore(1); // Skip one char before data
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Write the .pnm header.
            ///
            /// @param _type     Image type
            /// @param _rows     Nr of rows
            /// @param _cols     Nr of cols
            /// @param _maxval   Maxval
            /// @param comments Comments
            ////////////////////////////////////////////////////////////////////////////////////////////
            void write_header(const imtype _type, const arma::uword _rows, const arma::uword _cols, const int _maxval, const string comments)
            {
                type   = _type;
                rows   = _rows;
                cols   = _cols;
                maxval = _maxval;
                ofs << "P" << type << endl;
                ofs << "# " << comments << endl;
                ofs << cols << " " << rows << endl;
                if(type!=PBM_A || type!=PBM_B) ofs << maxval << endl;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Write the .pnm file.
            ///
            /// @returns true if success
            /// @param fname  File name
            /// @param _type  File name
            /// @param img    Image data
            /// @param info   File comments
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool write(std::string fname, const imtype _type, const arma::cube& img, const string info="")
            {
                // Open file
                ofs.open(fname.c_str(), ofstream::binary);
                if (!ofs.good())
                {
                    cout << "Could not open " << fname << endl;
                    return false;
                }
                else
                {
                    write_header(_type, img.n_rows, img.n_cols, img.max(),info);
                    //            get_info();

                    // Write data
                    if(type==PPM_A ) // Plain (ASCII )type
                    {
                        for(arma::uword r=0; r<rows; r++)
                        {
                            arma::uword i = 0;
                            for(arma::uword c=0; c<cols; c++)
                            {
                                ofs << img(r,c,0) << " " << img(r,c,1) << " " << img(r,c,2) << " "; // R G B
                                if(++i%5==0) ofs << endl; // Max len is 70 chars/line
                            }
                        }
                    }
                    else if(type==PPM_B)
                    {
                        for(arma::uword r=0; r<rows; r++)
                            for(arma::uword c=0; c<cols; c++)
                            {
                                unsigned char bb;
                                bb= static_cast<unsigned char>(img(r,c,0));   // R
                                ofs.write(reinterpret_cast<char*>(&bb),1);
                                bb= static_cast<unsigned char>(img(r,c,1));   // G
                                ofs.write(reinterpret_cast<char*>(&bb),1);
                                bb= static_cast<unsigned char>(img(r,c,2));   // B
                                ofs.write(reinterpret_cast<char*>(&bb),1);
                            }
                    }
                }
                ofs.close();
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Write the .pnm file.
            ///
            /// @returns true if success
            /// @param fname  File name
            /// @param _type  Image type
            /// @param img    Image data
            /// @param info   File comments
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool write(std::string fname, const imtype _type, arma::mat& img, const string info="")
            {
                // Open file
                ofs.open(fname.c_str(), ofstream::binary);
                if (!ofs.good())
                {
                    cout << "Could not open " << fname << endl;
                    return false;
                }
                else
                {
                    write_header(_type, img.n_rows, img.n_cols,img.max(),info);
                    //            get_info();

                    // Write data
                    if(type==PBM_A || type ==PGM_A ) // Plain (ASCII )type
                    {
                        arma::uword i=0;
                        for(arma::mat::iterator ii=img.begin(); ii!=img.end(); ++ii)
                        {
                            ofs << *ii << " ";
                            if(++i%11==0) ofs << endl; // Max len is 70 chars/line
                        }
                    }
                    else if(type == PBM_B) // Raw .pbm
                    {
                        bitset<8> b;
                        for(arma::uword r=0; r<rows; r++)
                            for(arma::uword c=0; c<cols; c++)
                            {
                                arma::uword ix = 7-(c%8);
                                b[ix] =  (img(r,c)>0);
                                if(ix==0 || c==cols-1)
                                {
                                    ofs.write(reinterpret_cast<char*>(&b),1);
                                    b.reset();
                                }
                            }
                    }
                    else if(type == PGM_B) // Raw .pgm
                    {
                        if(maxval<=255)
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned char bb= static_cast<unsigned char>(img(r,c));
                                    ofs.write(reinterpret_cast<char*>(&bb),1);
                                }
                        }
                        else
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned int bb;
                                    bb = ((static_cast<unsigned int>(img(r,c)))>> 8) & 0x00ff;   // Write MSB first
                                    ofs.write(reinterpret_cast<char*>(&bb),1);
                                    bb = static_cast<unsigned int>(img(r,c)) & 0x00ff;
                                    ofs.write(reinterpret_cast<char*>(&bb),1);
                                }
                        }
                    }

                }

                ofs.close();
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Prints header info
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////
            void get_info()
            {
                cout << "Type:  P"  << type   << endl;
                cout << "cols:   "  << cols   << endl;
                cout << "rows:   "  << rows   << endl;
                if(type==PGM_A || type==PGM_B) cout << "Maxval: "  << maxval << endl;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get nr of rows
            ///
            /// @returns number of rows
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::uword get_rows()
            {
                return rows;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get nr of cols
            ///
            /// @returns number of columns
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::uword get_cols()
            {
                return cols;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get maxval
            ///
            /// @returns Maximum value in image
            ////////////////////////////////////////////////////////////////////////////////////////////
            int get_maxval()
            {
                return maxval;
            }


            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Read image
            ///
            /// @returns true if success
            /// @param fname File name
            /// @param img   Image data
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool read(std::string fname, arma::cube& img)
            {
                // Open file
                ifs.open(fname.c_str(), ifstream::binary);
                if (!ifs.good())
                {
                    cout << "Could not open " << fname << endl;
                    return false;
                }
                else
                {
                    read_header();
                    //            get_info();

                    img.set_size(rows,cols,3);
                    arma::uword r = 0, c = 0;
                    // Get the data
                    if (type==PPM_A )  // Plain .PPM
                    {
                        string str;
                        arma::uword i=0;
                        while (ifs >> str && r<rows ) // Read until eof
                        {
                            // Remove comments
                            size_t mark = str.find_first_of("#");
                            if ( mark!= string::npos)
                            {
                                ifs.ignore(256, '\n');
                                str.erase(mark, 1);

                                if (str.empty()) continue;
                            }
                            int pix= atoi(str.c_str());  // Convert to int

                            img(r, c, i%3) = pix;
                            i++;
                            if(i%3==0)
                                if (++c == cols)
                                {
                                    c = 0;
                                    r++;
                                }
                        }
                    }
                    else if (type==PPM_B )  // Raw .PPM
                    {
                        for(arma::uword r=0; r<rows; r++)
                            for(arma::uword c=0; c<cols; c++)
                            {
                                unsigned char bb;
                                ifs.read(reinterpret_cast<char*>(&bb),1);    // R
                                img(r,c,0) = bb;
                                ifs.read(reinterpret_cast<char*>(&bb),1);    // G
                                img(r,c,1) = bb;
                                ifs.read(reinterpret_cast<char*>(&bb),1);    // B
                                img(r,c,2) = bb;
                            }
                    }
                }
                ifs.close();
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Read image
            ///
            /// @returns true if success
            /// @param fname File name
            /// @param img Image data
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool read(std::string fname, arma::mat& img)
            {
                // Open file
                ifs.open(fname.c_str(), ifstream::binary);
                if (!ifs.good())
                {
                    cout << "Could not open " << fname << endl;
                    return false;
                }
                else
                {
                    read_header();
                    get_info();

                    img.set_size(rows,cols);
                    arma::uword r = 0, c = 0;
                    // Get the data
                    if (type==PBM_A || type == PGM_A)  // Plain .PGM or .PBM
                    {
                        string str;
                        while (ifs >> str && r<rows ) // Read until eof
                        {
                            // Remove comments
                            size_t mark = str.find_first_of("#");
                            if ( mark!= string::npos)
                            {
                                ifs.ignore(256, '\n');
                                str.erase(mark, 1);

                                if (str.empty()) continue;
                            }
                            int pix= atoi(str.c_str());  // Convert to int
                            img(r, c) = pix;

                            if (++c == cols)
                            {
                                c = 0;
                                r++;
                            }
                        }
                    }
                    else if(type== PBM_B) // Raw PBM
                    {
                        unsigned char ch;
                        while (ifs.read(reinterpret_cast<char*>(&ch),1) && r<rows)  // Read until eof
                        {
                            bitset<8> pix(ch);
                            for(int b=7; b>=0; b--)
                            {
                                img(r,c) = pix[b];
                                if (++c >= cols)
                                {
                                    c = 0;
                                    r++;
                                    break;
                                }
                            }
                        }
                    }
                    else if(type==PGM_B) // Raw PGM
                    {
                        if(maxval<=255)
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned char bb;
                                    ifs.read(reinterpret_cast<char*>(&bb),1);
                                    img(r,c) = bb;
                                }
                        }
                        else
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned char bb[2];
                                    ifs.read(reinterpret_cast<char*>(bb),2);
                                    img(r,c) = (bb[0]<<8)+bb[1];
                                }
                        }
                    }

                }
                ifs.close();
                return true;
            }

    };

    /// @}
}
#endif
