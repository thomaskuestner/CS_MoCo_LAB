/*
 * Copyright (c) 2010-2012, Andreas Kotowicz
 *
 *
 * ideas for this code are from:
 * http://download.oracle.com/javase/tutorial/networking/sockets/clientServer.html
 * http://www.mathworks.com/matlabcentral/fileexchange/24594-parfor-progress-monitor
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */ 

import java.io.*;
import java.net.*;

/* 
 * THIS CLASS WILL only BE USED IF
 *
 *
 *  1) you call it from the command line
 *
 * OR
 *
 *  2) you have 'matlabpool' running
 *
 *
 */


public class ParforProgressClient2 {
    
    // server hostname and port number
    private String sHost;    
    private int sPort;
    
    // show debugging information?
    private int DEBUG = 0;
    
    // timeout variables
    private final int TIME_TO_SLEEP_AFTER_TIMEOUT = 250;
    private int TimeoutCounter = 0;
    private int TimeoutCounterWarningThreshold = 100;
    
    // error indicators
    private boolean TimeoutWarningShown = false;
    private boolean error = false;
    private boolean UnknownHostErrorShown = false;
    private boolean SocketConnectionOpenErrorShown = false;
    private boolean SocketConnectionCloseErrorShown = false;
    
    
    /**
     * Create a "client" which runs on the remote lab and connects to our 
     * server via the network
     */
    public static ParforProgressClient2 createClient( String host, int port, int do_debug ) throws IOException {
        ParforProgressClient2 myclient = new ParforProgressClient2();
        myclient.setup(host, port, do_debug);
        return myclient;
    }
    
    /**
     * setup the parameters
     */
    public void setup( String hostname, int portnumber, int do_debug) {
        sPort = portnumber;
        sHost = hostname;
        DEBUG = do_debug;
        
        if (DEBUG > 1)
            System.out.println("ParforProgressClient2: setup() complete.");
    }
    
    
    /* 
     * this method will be used if matlabpool is ON
     */
    public void increment() {

        if (DEBUG > 1)
            System.out.println("ParforProgressClient2: calling increment().");
        
        Socket clientSocket = null;
        
        // try to connect to server.
        try {
            
            // open new socket
            clientSocket = new Socket(sHost, sPort);
            
            // maybe this helps with ultra short connection times ...
            clientSocket.setReuseAddress(true);
            
            // takes care of connections being left in 'TIME_WAIT' state
            clientSocket.setSoLinger(true, 0);
            
        } catch (UnknownHostException e) {
            error = true;
            if (DEBUG > 0 && UnknownHostErrorShown == false) {
                System.err.println("ParforProgressClient2: Host not found: " + sHost);
                System.err.println(e);                
                UnknownHostErrorShown = true;
            }
        } catch (IOException e) {
            error = true;
            if (DEBUG > 0 && SocketConnectionOpenErrorShown == false) {
                System.err.println("ParforProgressClient2: Couldn't setup socket connection to: " + sHost);
                System.err.println(e);
                SocketConnectionOpenErrorShown = true;
            }
        }
        
        // pause execution in case of previous connection error.
        if (error == true) {
            try {
                Thread.sleep(TIME_TO_SLEEP_AFTER_TIMEOUT);
                error = false;
                // update TimeoutCounter only up to threshold - otherwise
                // it might explode.
                if (TimeoutCounter <= TimeoutCounterWarningThreshold) {
                    TimeoutCounter++;
                } else {
                    if (TimeoutWarningShown == false) {
                        System.err.println("ParforProgressClient2: Your single loop iteration is too easy to compute. Try to update the progress only each 10 loops or so.");
                        TimeoutWarningShown = true;
                    }
                }
            }
            catch (InterruptedException f) {
            }
        }
        
        // close connection if possible.
        if (clientSocket != null) {
            try {
                clientSocket.close();
            } catch (IOException e) {
                if (DEBUG > 0 && SocketConnectionCloseErrorShown == false) {
                    System.err.println("ParforProgressClient2: Problem closing socket.");
                    System.err.println(e);
                    SocketConnectionCloseErrorShown = true;
                }
            }
        }

        /* don't do 
         *
         *      System.exit(1);
         *
         * here, otherwise the matlab parfor loop will break!
         */
        
    }
    
    
    /**
     * Nothing for us to do here. done() only exists for compatibility with
     * ParforProgressServer2 (we never know whether the user is 
     * running with or without threads).
     */
    public void done() {
        if (TimeoutCounter > 0)
            System.err.println("ParforProgressClient2: timeouts: " + TimeoutCounter);
        if (DEBUG > 1)
            System.err.println("ParforProgressClient2: done() called.");
    }    
    
    
    public static void main(String[] args) {

        if (args.length < 2) {
            System.err.println("ParforProgressClient2: Please provide at least 2 input arguments (hostname, port).");
            System.exit(1);
        }
        
        String hostname = args[0];
        int port = 10;
        int do_debug = 0;
        
        if (args.length > 1) {
            try {
                port = Integer.parseInt(args[1]);
            } catch (NumberFormatException e) {
                System.err.println("ParforProgressClient2: 2nd argument must be an integer.");
                System.exit(1);
            }
        }
        
        if (args.length > 2) {
            try {
                do_debug = Integer.parseInt(args[2]);
            } catch (NumberFormatException e) {
                System.err.println("ParforProgressClient2: 3rd argument must be an integer.");
                System.exit(1);
            }
        }
        
        ParforProgressClient2 myclient = new ParforProgressClient2();
        myclient.setup(hostname, port, do_debug);
        myclient.increment();
    }
    
}
