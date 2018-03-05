% load xml file with xml toolbox from Marc Molinari (MATLAB file exchange)


% quick 'n dirty - don't change order in XML file!!


function sTrafo = fLoadXML(sFilename, sMeasPara)
    sXML = xml_parseany(fileread(sFilename));
    
    lambda                  = str2num(sXML.gadget{1,2}.property{1,1}.value{1,1}.CONTENT);
    lambdaESPReSSo          = str2num(sXML.gadget{1,2}.property{1,2}.value{1,1}.CONTENT);
    algorithm               = str2num(sXML.gadget{1,2}.property{1,3}.value{1,1}.CONTENT);
    iNOUTER                 = str2num(sXML.gadget{1,2}.property{1,4}.value{1,1}.CONTENT);
    iNINNER                 = str2num(sXML.gadget{1,2}.property{1,5}.value{1,1}.CONTENT);
    Kernel_FFT_dim          = str2num(sXML.gadget{1,2}.property{1,6}.value{1,1}.CONTENT);
    FFT_Sparse              = str2num(sXML.gadget{1,2}.property{1,7}.value{1,1}.CONTENT);
    DCT_Sparse              = str2num(sXML.gadget{1,2}.property{1,8}.value{1,1}.CONTENT);
    PCA_Sparse              = str2num(sXML.gadget{1,2}.property{1,9}.value{1,1}.CONTENT);
    kSpaceTrafo             = str2num(sXML.gadget{1,2}.property{1,10}.value{1,1}.CONTENT);
    Transform_fftBA_dim     = str2num(sXML.gadget{1,2}.property{1,11}.value{1,1}.CONTENT);
    kSpaceOut               = str2num(sXML.gadget{1,2}.property{1,12}.value{1,1}.CONTENT);
    CS_ESPReSSo             = str2num(sXML.gadget{1,2}.property{1,13}.value{1,1}.CONTENT);
    
    sTrafo = struct('lambda', lambda, 'lambdaESPReSSo', lambdaESPReSSo,...
                    'algorithm', algorithm, 'iNOUTER', iNOUTER, 'iNINNER', iNINNER, 'Kernel_FFT_dim', Kernel_FFT_dim,...
                    'FFT_Sparse', FFT_Sparse, 'DCT_Sparse', DCT_Sparse,'PCA_Sparse',PCA_Sparse,...
                    'kSpaceTrafo', kSpaceTrafo,'Transform_fftBA_dim', Transform_fftBA_dim,'kSpaceOut', kSpaceOut,...
                    'CS_ESPR',CS_ESPReSSo);
                
    sTrafo.CSFullySampled       = sMeasPara.CSFullySampled;
    sTrafo.ESPReSSoDirection    = sMeasPara.ESPReSSoDirection;
    sTrafo.dataset              = 2;
    sTrafo.dPartialFourierVal   = sMeasPara.ESPReSSo;
end