function hF1 = fEvalRegistrationGray(dIRef, dIMove, dFx, dFy)
    
    h.iQUIVERFACTOR = 6;
    dSCALEFACTOR = 1;

    dIRef  = double(dIRef );
    dIMove = double(dIMove);
    dIRef  = dIRef  - min(dIRef (:));
    dIMove = dIMove - min(dIMove(:));
    dIRef  = dIRef ./max(dIRef (:)) .* dSCALEFACTOR;
    dIMove = dIMove./max(dIMove(:)) .* dSCALEFACTOR;
    dIRef (dIRef  > 1) = 1;
    dIMove(dIMove > 1) = 1;
    
    h.dIZero = zeros(size(dIRef, 1), size(dIRef, 2));
    oldMax = max(dIRef(:));
%     dIRef = dIRef .* dSCALEFACTOR ;

    hF1 = figure;
    h.hI = imshow(dIRef(:,:,1),[0 oldMax]);
%     h.hI = imshow(cat(3, dIRef(:,:,1), dIMove(:,:,1), h.dIZero));
    h.hA = gca;
    set(h.hA, 'Position', [0 0 1 1]);
    set(hF1, 'Position', [0 0 size(dIMove, 2).*4 size(dIMove, 1).*4]);
    set(hF1, 'WindowScrollWheelFcn', @fScroll);
    movegui('center');
    hold on
    
    h.dIRef = dIRef;
    h.dIMove = dIMove;
    h.dFx = dFx;
    h.dFy = dFy;
%    h.dFz = dFz;
    h.iActive = 1;
    
    h.iX = 1:h.iQUIVERFACTOR:size(h.dFy, 2);
    h.iY = 1:h.iQUIVERFACTOR:size(h.dFy, 1);
        
    h.hQ = quiver(h.iX, h.iY, h.dFx(1:h.iQUIVERFACTOR:end, 1:h.iQUIVERFACTOR:end, 1), ...
                              h.dFy(1:h.iQUIVERFACTOR:end, 1:h.iQUIVERFACTOR:end, 1), 0);
    
    set(h.hQ, 'Linewidth', 1.5, 'Color', 'y');
    guidata(hF1, h);
end

function fScroll(hObject, eventdata, handles)

    h = guidata(hObject);
    if eventdata.VerticalScrollCount < 0
        h.iActive = max([1 h.iActive - 1]);
    else
        h.iActive = min([size(h.dIRef, 3) h.iActive + 1]);
    end
%     set(h.hI, 'CData', cat(3, h.dIMove(:,:,h.iActive), h.dIRef(:,:,h.iActive), h.dIZero));
    set(h.hI, 'CData', cat(3, h.dIRef(:,:,h.iActive), h.dIMove(:,:,h.iActive), h.dIZero));
%     set(h.hI, 'CData', h.dIRef(:,:,h.iActive));
    set(h.hQ, 'UData', h.dFx(1:h.iQUIVERFACTOR:end, 1:h.iQUIVERFACTOR:end, h.iActive), ...
              'VData', h.dFy(1:h.iQUIVERFACTOR:end, 1:h.iQUIVERFACTOR:end, h.iActive));
    guidata(hObject, h);
end