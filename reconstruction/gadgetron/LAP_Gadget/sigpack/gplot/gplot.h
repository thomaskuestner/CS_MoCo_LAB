// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_GPLOT_H
#define SP_GPLOT_H
namespace sp
{
    ///
    /// @defgroup gplot GPlot
    /// \brief Collection of Gnuplot functions
    /// @{

    ///
    /// \brief Gnuplot class.
    ///
    /// Implements a class for streaming data to Gnuplot using a pipe.
    /// Inspiration from https://code.google.com/p/gnuplot-cpp/
    ///
    /// Verified with Gnuplot 4.6.5 for Win64 and Linux.
    /// \note In Windows only one class is allowed. Using multiple figures are controlled by a figure number. In Linux we may use one instance per figure.
    ///
    class gplot
    {
        private:
            FILE*           gnucmd;          ///< File handle to pipe
            std::string     linestyle;       ///< Linestyle string as in Gnuplot doc.
            std::string     term;

        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ///
            /// Opens a pipe to gnuplot program. Make sure it is found/accessable by the system.
            ////////////////////////////////////////////////////////////////////////////////////////////
            gplot()
            {
#if defined(WIN32)
                gnucmd = _popen("gnuplot -persist 2> NUL","w");
                term = "win";
#elif defined(unix)
                //gnucmd = popen("gnuplot -persist &> /dev/null","w");
                gnucmd = popen("gnuplot -persist","w");
                term = "x11";
                //#elif defined(_APPLE_)
                //            gnucmd = popen("gnuplot -persist &> /dev/null","w");
                //#define term "aqua"
#else
#error Only Windows and Linux/Unix is supported so far!
#endif
                if(!gnucmd)
                {
                    err_handler("Could not start gnuplot");
                }
                setvbuf(gnucmd, NULL, _IOLBF, 512);

                // Set global params
                linestyle = "lines";
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~gplot()
            {
#if defined(WIN32)
                _pclose(gnucmd);
#elif defined(unix)
                pclose(gnucmd);
#endif
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Send command to Gnuplot pipe.
            /// @param cmdstr  Command string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void send2gp(const char* cmdstr)
            {
                std::string s_in(cmdstr);
                std::string tmp=s_in+"\n";
                std::fputs(tmp.c_str(), gnucmd );
                //        std::cout << tmp.c_str() << std::endl;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Sets the active figure.
            /// @param fig  Figure number
            ////////////////////////////////////////////////////////////////////////////////////////////
            void figure(const int fig)
            {
                std::ostringstream tmp_s;
                tmp_s << "set term " << term << " " << fig;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                send2gp("reset");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Configure the figure used Windows environment.
            /// @param fig     Figure number
            /// @param name    Window name
            /// @param x       x position of upper left corner
            /// @param y       y position of upper left corner
            /// @param width   width of window
            /// @param height  height of window
            ////////////////////////////////////////////////////////////////////////////////////////////
            void window(const int fig, const char* name,const int x,const int y,const int width,const int height)
            {
                std::ostringstream tmp_s;
                tmp_s << "set term " << term << " " << fig << " title \"" << name << "\" position " << x << "," << y << " size " << width << "," << height;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                send2gp("reset");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Configure the figure/window - used in Linux environment where no figure numbers are needed.
            /// @param name    Window name
            /// @param x       x position of upper left corner
            /// @param y       y position of upper left corner
            /// @param width   width of window
            /// @param height  height of window
            ////////////////////////////////////////////////////////////////////////////////////////////
            void window(const char* name,const int x,const int y,const int width,const int height)
            {
                window(0,name,x,y,width,height);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Close window
            ////////////////////////////////////////////////////////////////////////////////////////////
            void close_window(void)
            {
                std::ostringstream tmp_s;
                tmp_s << "set term " << term << " close";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                send2gp("reset");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set linestyle
            /// @param style  lines | points | linespoints | dots | steps
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_linestyle(const char* style)
            {
                linestyle = style;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set label for X-axis.
            /// @param label label string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void xlabel(const char* label)
            {
                std::ostringstream tmp_s;
                tmp_s << "set xlabel \"" << label << "\" ";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set label for X-axis.
            /// @param label label string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ylabel(const char* label)
            {
                std::ostringstream tmp_s;
                tmp_s << "set ylabel \"" << label << "\" ";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set label at position x,y.
            /// @param x x value
            /// @param y y value
            /// @param label label string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void label(const double& x, const double& y, const char* label)
            {
                std::ostringstream tmp_s;
                tmp_s << "set label \"" << label << "\" at " << x << "," << y;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set windowtitle.
            /// @param name title string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void title(const char* name)
            {
                std::ostringstream tmp_s;
                tmp_s << "set title \"" << name << " \" ";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set X-axis range.
            /// @param xmin xmin
            /// @param xmax xmax
            ////////////////////////////////////////////////////////////////////////////////////////////
            void xlim(const double xmin, const double xmax)
            {
                std::ostringstream tmp_s;
                tmp_s << "set xrange [" << xmin << ":" << xmax << "]";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set Y-axis range.
            /// @param ymin ymin
            /// @param ymax ymax
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ylim(const double ymin, const double ymax)
            {
                std::ostringstream tmp_s;
                tmp_s << "set yrange [" << ymin << ":" << ymax << "]";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot y vs. x.
            /// @param x x vector
            /// @param y y vector
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot_str2( const arma::vec& x, const arma::vec& y)
            {
                std::ostringstream tmp_s;
                for(arma::uword n=0; n<y.size(); n++)
                {
                    tmp_s << x(n) << " " << y(n);
                    std::string s = tmp_s.str();
                    send2gp(s.c_str());
                    tmp_s.str(""); // Clear buffer
                }
                send2gp("e");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot multiple y vs. x
            /// @param x      x vector
            /// @param y      y matrix with data in rows
            /// @param label  data label prefix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot( const arma::vec& x, const arma::mat& y, const char* label="")
            {
                arma::uword R =  y.n_rows;
                std::string label_s(label);
                std::ostringstream tmp_s;
                send2gp("set key noautotitle");
                send2gp("set grid");
                if (label_s.empty())
                {
                    tmp_s << "plot '-' with " << linestyle;
                    for(arma::uword r=1; r<R; r++)
                        tmp_s << " ,'-' with " << linestyle;
                }
                else
                {
                    tmp_s << "plot '-' title \"" << label_s << "1\" with " << linestyle;
                    for(arma::uword r=1; r<R; r++)
                        tmp_s << " ,'-' title \"" << label_s << r+1 << "\" with " << linestyle;
                }
                std::string s = tmp_s.str();
                send2gp(s.c_str());

                for(arma::uword r=0; r<y.n_rows; r++)
                {
                    arma::vec line = arma::conv_to<arma::vec>::from(y.row(r));
                    plot_str2(x,line);
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot multiple y
            /// Each line plot is a row in matrix
            /// @param y     y matrix with data in rows
            /// @param label data label prefix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot( const arma::mat& y, const char* label="")
            {
                arma::vec t = arma::linspace(1,y.n_cols,y.n_cols);
                plot(t,y,label);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot y vs. x. with label
            /// @param x     x vector
            /// @param y     y vector
            /// @param label plot label
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot( const arma::vec& x, const arma::vec& y, const char* label="")
            {
                std::ostringstream tmp_s;
                std::string label_s(label);

                if (label_s.empty())
                {
                    send2gp("set nokey");
                    tmp_s << "plot '-' with " << linestyle;
                }
                else
                {
                    send2gp("set key noautotitle");
                    tmp_s << "plot '-' title \"" << label_s << "\" with " << linestyle;
                }
                send2gp("set grid");
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                plot_str2(x,y);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot y with label
            /// @param y     y vector
            /// @param label plot label
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot(const arma::vec& y, const char* label="")
            {
                arma::vec t(y.size());
                t = arma::linspace(1,y.size(),y.size());
                plot(t,y,label);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot y1 and y2 vs. x with label
            /// @param x      x vector
            /// @param y1     y1 vector
            /// @param y2     y2 vector
            /// @param label1 plot label 1
            /// @param label2 plot label 2
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot( const arma::vec& x, const arma::vec& y1, const arma::vec& y2, const char* label1, const char* label2)
            {
                std::ostringstream tmp_s;
                send2gp("set key noautotitle");
                send2gp("set grid");
                tmp_s << "plot '-' title \"" << label1 << "\" with " << linestyle << ", '-' title \"" << label2 << "\" with " << linestyle;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                plot_str2(x,y1);
                plot_str2(x,y2);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot y1 and y2 with label
            /// @param y1     y1 vector
            /// @param y2     y2 vector
            /// @param label1 plot label 1
            /// @param label2 plot label 2
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot(const arma::vec& y1, const arma::vec& y2, const char* label1, const char* label2)
            {
                arma::vec t(y1.size());
                t = arma::linspace(1,y1.size(),y1.size());
                plot(t,y1,y2,label1,label2);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Scatter plot
            /// @param x    x vector
            /// @param y    y vector
            /// @param label plot label
            ////////////////////////////////////////////////////////////////////////////////////////////
            void scatter( const arma::vec& x, const arma::vec& y, const char* label="")
            {
                set_linestyle("points");
                plot(x,y,label);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Scatter plot - dual
            /// @param x1     x1 vector
            /// @param y1     y1 vector
            /// @param label1 plot label1
            /// @param x2     x2 vector
            /// @param y2     y2 vector
            /// @param label2 plot label2
            ////////////////////////////////////////////////////////////////////////////////////////////
            void scatter( const arma::vec& x1, const arma::vec& y1, const char* label1,
                          const arma::vec& x2, const arma::vec& y2, const char* label2)
            {
                set_linestyle("points");
                std::ostringstream tmp_s;
                send2gp("set key noautotitle");
                send2gp("set grid");
                tmp_s << "plot '-' title \"" << label1 << "\" with " << linestyle << ", '-' title \"" << label2 << "\" with " << linestyle;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                plot_str2(x1,y1);
                plot_str2(x2,y2);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot mat as image
            /// @param x     x matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void image( const arma::mat& x)
            {
                std::ostringstream tmp_s;
                xlim(-0.5,x.n_cols-0.5);
                ylim(x.n_rows-0.5,-0.5);
                tmp_s.str(""); // Clear buffer
                tmp_s << "plot '-' matrix with image";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                for(arma::uword r=0; r<x.n_rows; r++)
                {
                    tmp_s.str("");  // Clear buffer
                    for(arma::uword c=0; c<x.n_cols; c++)
                    {
                        tmp_s << x(r,c) << " " ;
                    }
                    s = tmp_s.str();
                    send2gp(s.c_str());
                }
                send2gp("e");
                send2gp("e");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot cube as image
            /// @param x     x matrix (R,G,B)
            ////////////////////////////////////////////////////////////////////////////////////////////
            void image( const arma::cube& x)
            {
                std::ostringstream tmp_s;
                xlim(-0.5,x.n_cols-0.5);
                ylim(x.n_rows-0.5,-0.5);
                tmp_s.str(""); // Clear buffer
                tmp_s << "plot '-' with rgbimage";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                for(arma::uword r=0; r<x.n_rows; r++)
                {
                    tmp_s.str("");  // Clear buffer
                    for(arma::uword c=0; c<x.n_cols; c++)
                    {
                        tmp_s << " " << c << " " << r << " " << x(r,c,0) << " " << x(r,c,1) << " " << x(r,c,2) << std::endl;
                    }
                    std::string s = tmp_s.str();
                    send2gp(s.c_str());
                }
                send2gp("e");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot mat as mesh
            /// @param x     x matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void mesh( const arma::mat& x)
            {
                std::ostringstream tmp_s;
                send2gp("unset key");
                send2gp("set hidden3d");
                tmp_s << "splot '-' with lines";
                std::string s = tmp_s.str();
                send2gp(s.c_str());

                for(arma::uword r=0; r<x.n_rows; r++)
                {
                    for(arma::uword c=0; c<x.n_cols; c++)
                    {
                        tmp_s.str("");  // Clear buffer
                        tmp_s << r << " " << c << " "<< x(r,c);
                        std::string s = tmp_s.str();
                        send2gp(s.c_str());
                    }
                    send2gp("");
                }
                send2gp("e");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot mat as surf
            /// @param x     x matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void surf( const arma::mat& x)
            {
                send2gp("set pm3d");
                mesh(x);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set linestyle to Matlab 'parula' NB! doesn't work with X11 -terminal
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_parula_line(void)
            {
                send2gp("set style line  1 lt 1 lc rgb '#352a87' "); // blue
                send2gp("set style line  2 lt 1 lc rgb '#0f5cdd' "); // blue
                send2gp("set style line  3 lt 1 lc rgb '#1481d6' "); // blue
                send2gp("set style line  4 lt 1 lc rgb '#06a4ca' "); // cyan
                send2gp("set style line  5 lt 1 lc rgb '#2eb7a4' "); // green
                send2gp("set style line  6 lt 1 lc rgb '#87bf77' "); // green
                send2gp("set style line  7 lt 1 lc rgb '#d1bb59' "); // orange
                send2gp("set style line  8 lt 1 lc rgb '#fec832' "); // orange
                send2gp("set style line  9 lt 1 lc rgb '#f9fb0e' "); // yellow
                send2gp("set style line 11 lt 1 lc rgb '#0072bd' "); // blue
                send2gp("set style line 12 lt 1 lc rgb '#d95319' "); // orange
                send2gp("set style line 13 lt 1 lc rgb '#edb120' "); // yellow
                send2gp("set style line 14 lt 1 lc rgb '#7e2f8e' "); // purple
                send2gp("set style line 15 lt 1 lc rgb '#77ac30' "); // green
                send2gp("set style line 16 lt 1 lc rgb '#4dbeee' "); // light-blue
                send2gp("set style line 17 lt 1 lc rgb '#a2142f' "); // red
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to Matlab 'jet'
            ///  See https://github.com/Gnuplotting/gnuplot-palettes
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_jet_palette(void)
            {
                char str[] ="set palette defined ( \
              0 '#000090',\
              1 '#000fff',\
              2 '#0090ff',\
              3 '#0fffee',\
              4 '#90ff70',\
              5 '#ffee00',\
              6 '#ff7000',\
              7 '#ee0000',\
              8 '#7f0000')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to Matlab 'parula'
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_parula_palette(void)
            {
                char str[] ="set palette defined (\
              0 '#352a87',\
              1 '#0363e1',\
              2 '#1485d4',\
              3 '#06a7c6',\
              4 '#38b99e',\
              5 '#92bf73',\
              6 '#d9ba56',\
              7 '#fcce2e',\
              8 '#f9fb0e')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to 'cool-warm'
            // See http://www.kennethmoreland.com/color-advice/
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_coolwarm_palette(void)
            {
                char str[] = "set palette defined (\
              0 '#5548C1', \
              1 '#7D87EF', \
              2 '#A6B9FF', \
              3 '#CDD7F0', \
              4 '#EBD1C2', \
              5 '#F3A889', \
              6 '#DE6A53', \
              7 '#B10127')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to 'black body'
            // See http://www.kennethmoreland.com/color-advice/
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_blackbody_palette(void)
            {
                char str[] = "set palette defined (\
              0 '#000000', \
              1 '#2B0F6B', \
              2 '#5D00CB', \
              3 '#C60074', \
              4 '#EB533C', \
              5 '#F59730', \
              6 '#E9D839', \
              7 '#FFFFFF')";
                send2gp(str);
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Save plot to file.
            /// @param name filename
            ///
            /// Extensions that are supported:
            /// - png
            /// - ps
            /// - eps
            /// - tex
            /// - pdf
            /// - svg
            /// - emf
            /// - gif
            ///
            /// \note When 'latex' output is used the '\' must be escaped by '\\\\' e.g set_xlabel("Frequency $\\\\omega = 2 \\\\pi f$")
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_output(const char* name)
            {
                std::string name_s(name);
                size_t found = name_s.find_last_of(".");
                std::string ext;
                ext = name_s.substr(found + 1);
                std::ostringstream tmp_s;

                if (ext.compare("png")==0)
                {
                    tmp_s << "set terminal pngcairo enhanced font 'Verdana,10'";
                }
                else if (ext.compare("ps") == 0)
                {
                    tmp_s << "set terminal postscript enhanced color";
                }
                else if (ext.compare("eps") == 0)
                {
                    tmp_s << "set terminal postscript eps enhanced color";
                }
                else if (ext.compare("tex") == 0)
                {
                    tmp_s << "set terminal cairolatex eps color enhanced";
                }
                else if (ext.compare("pdf") == 0)
                {
                    tmp_s << "set terminal pdfcairo color enhanced";
                }
                else if (ext.compare("svg") == 0)
                {
                    tmp_s << "set terminal svg enhanced";
                }
                else if (ext.compare("emf") == 0)
                {
                    tmp_s << "set terminal emf color enhanced";
                }
                else if (ext.compare("gif") == 0)
                {
                    tmp_s << "set terminal gif enhanced";
                }
                //else if (ext.compare("jpg") == 0)
                //{
                //	tmp_s << "set terminal jpeg ";
                //}
                else
                {
                    tmp_s << "set terminal " << term;
                }
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                tmp_s.str("");  // Clear buffer
                tmp_s << "set output '" << name_s << "'";
                s = tmp_s.str();
                send2gp(s.c_str());
            }


            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Restore output terminal.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void restore_output(void)
            {
                std::ostringstream tmp_s;
                tmp_s << "set terminal " << term;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set output terminal.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_term(const char* ttype)
            {
                std::ostringstream tmp_s;
                tmp_s << "set terminal " << ttype;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

    }; // End Gnuplot Class


    /// @}

} // end namespace
#endif
