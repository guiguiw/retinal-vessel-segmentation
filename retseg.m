function CLEANED = retseg(INNAME,MASKNAME)
%RETINAL SEGMENATION 
%   CLEANED = retseg(INNAME,MASKNAME) segments blood vessels from full
%   color retinal image INNAME. MASKNAME is the corresponding mask image
%   as provided in the DRIVE database. Provided code could be modified
%   slightly to create mask image from input image.
%
%   INNAME and MASKNAME are to be passed as strings. 
%
%   Output CLEANED is a logical matrix, bright vessels indicate pixels
%   found to be vessel pixels.
%
%   Runtime on school computers in ELW B326: ~ 25 mins
%   Runtime on Brandon's Laptop: ~40 mins
%   CREATED USING MATLAB 2014 DOCUMENTATION
%
% example to run code using provided images
% CLEANED = retseg('21_training.tif','21_training_mask.gif');
% CLEANED = retseg('01_test.tif','01_test_mask.gif');
% CLEANED = retseg('04_test.tif','04_test_mask.gif');
% Link to DRIVE database: http://www.isi.uu.nl/Research/Databases/DRIVE/

starter = tic;          %Calculates run time

%thresholds from paper
Pphi = 135;
Pd = 9;
Psize = 37;
Pfac = 2.0;
Pcontrast = 1.07;
Ts = 60:2:220; %PAPER USED STEP SIZE OF 2 

%inputs - to be removed and replaced with automated code
% inname = '21_training.tif';
maskimage = imread(MASKNAME);
inputimage = imread(INNAME);            %reads the input image
greenimage = inputimage(:,:,2);         %pulls out the green layer only
[rowmax,colmax] = size(greenimage);

% variable intialization
contrastimage = zeros(rowmax,colmax,length(Ts));
threshmap = 255*ones(rowmax,colmax);
SEGMENTED = zeros(size(greenimage));
masklogical=maskimage>1;
disp('Image name');
disp(INNAME)
%The algorithm thresholds an image multiple times, this for loop iterates
%through each threshold value.

LTs = sprintf('%d Layers Total', length(Ts));
disp(LTs)
disp('');
disp('Current Layer');
for layer = 1:length(Ts)
    
    disp(layer)
        
    %------------------------Thresholding section--------------------------
    %This is used to find pixels who's value is greater than or less than a
    %certain luminance threshold. Output threshmap is binary
    %could be a function with inputs (workingimage,Ts) and would return
    %threshmap, a binary matrix/image
            
    for rowdex = 1:rowmax
        for coldex = 1:colmax
            if greenimage(rowdex,coldex) < Ts(layer) 
                    threshmap(rowdex,coldex) = 0;
            end %end if
        end     %end for (cols)
    end         %end for (rows)
    
    
    %----------------------Euclidean Distance section---------------------
    %could be a function with input (threshmap) and outputs
    %[distmat,rowcoord,colcoord]
    %distmat is a matrix who's value is the distance to the closest
    %background pixel
    %rowcoord and colcoords show the row and column index respectively of
    %said pixel
    [distmat, linindex] = bwdist(threshmap);
    
    %breaks linindex into meaningful row and column coordinates
    [rowcoord,colcoord] = ind2sub(size(linindex),linindex);
    
        
    %---------------------------Quick Filter-----------------------------
    %pseudocode 
    % if length(vector(p e_p)) > P_d / ( sqrt(2) * sqrt(1-cos(P_phi)) ) 
    % then you toss p. 
    % not implemented due to worse results, and since it only affected
    % computational efficiency, we deemed it okay to remove
    
    %turns boundary/background pixels into true 0 value
    distmasked = distmat.*masklogical;
    %which aids run time
    
    %---------------------------Skeleton code----------------------------
    %finds the skeleton outline and disregards anything that isn't part of
    %it
    
    for rowdex = 2:rowmax-1 %fortunately, we don't have to pad due to the meaningful data being within this border
        for coldex = 2:colmax-1
            if distmasked(rowdex,coldex) > 0      %if p is a candidate 

                %-----------------Vector intitialization-------------------

                %current location within layer
                P = [rowdex,coldex];
                
                %corresponding background pixel
                EP = [rowcoord(rowdex,coldex),colcoord(rowdex,coldex)];
                EP = double(EP);

                %neighbors (N8)
                %MATLAB vectors of spatial points. N8 neighbors of P closest
                %background pixels.

                EN1 = [rowcoord(rowdex-1,coldex-1),colcoord(rowdex-1,coldex-1)];    %top left
                EN2 = [rowcoord(rowdex,coldex-1),colcoord(rowdex,coldex-1)];        %top
                EN3 = [rowcoord(rowdex+1,coldex-1),colcoord(rowdex+1,coldex-1)];    %top right
                EN4 = [rowcoord(rowdex+1,coldex),colcoord(rowdex+1,coldex)];        %right
                EN5 = [rowcoord(rowdex+1,coldex+1),colcoord(rowdex+1,coldex+1)];    %bottom right
                EN6 = [rowcoord(rowdex,coldex+1),colcoord(rowdex,coldex+1)];        %bottom
                EN7 = [rowcoord(rowdex-1,coldex+1),colcoord(rowdex-1,coldex+1)];    %bottom left
                EN8 = [rowcoord(rowdex-1,coldex),colcoord(rowdex-1,coldex)];        %left
                
                %for type requirements later
                EN1 = double(EN1);
                EN2 = double(EN2);
                EN3 = double(EN3);
                EN4 = double(EN4);
                EN5 = double(EN5);
                EN6 = double(EN6);
                EN7 = double(EN7);
                EN8 = double(EN8);
                                               
                %vector addition
                PEP = EP-P;
                PEN1 = EN1-P;
                PEN2 = EN2-P;
                PEN3 = EN3-P;
                PEN4 = EN4-P;
                PEN5 = EN5-P;
                PEN6 = EN6-P;
                PEN7 = EN7-P;
                PEN8 = EN8-P;

                %----------------------Angle Pruning-------------------------
                %can be turned into a function with
                %(PEP,PEN1,PEN2,PEN3,PEN4,PEN5,PEN6,PEN7,PEN8)
                %as arguments.
                %output would be the max angle
                %pseudocode                
                %for each neighbor, 
                %phi = max((180/pi)*acos((vector(p e_p) dot vector(p e_n))/norm(vector(p e_p) * norm(vector(p e_n)))

                phi(1) = 180/pi * acos((dot(PEP,PEN1))/norm(PEP)*norm(PEN1));
                phi(2) = 180/pi * acos((dot(PEP,PEN2))/norm(PEP)*norm(PEN2));
                phi(3) = 180/pi * acos((dot(PEP,PEN3))/norm(PEP)*norm(PEN3));
                phi(4) = 180/pi * acos((dot(PEP,PEN4))/norm(PEP)*norm(PEN4));
                phi(5) = 180/pi * acos((dot(PEP,PEN5))/norm(PEP)*norm(PEN5));
                phi(6) = 180/pi * acos((dot(PEP,PEN6))/norm(PEP)*norm(PEN6));
                phi(7) = 180/pi * acos((dot(PEP,PEN7))/norm(PEP)*norm(PEN7));
                phi(8) = 180/pi * acos((dot(PEP,PEN8))/norm(PEP)*norm(PEN8));

                phimax = max(phi);

                if phimax > Pphi          %if phi is wide enough
                    %-------------------Distance Pruning-------------------
                    %can be turned into a funtion with
                    %(PEP,PEN1,PEN2,PEN3,PEN4,PEN5,PEN6,PEN7,PEN8)
                    %as arguments
                    
                    %output would be the maximum distance

                    dist(1) = norm(PEP-PEN1);
                    dist(2) = norm(PEP-PEN2);
                    dist(3) = norm(PEP-PEN3);
                    dist(4) = norm(PEP-PEN4);
                    dist(5) = norm(PEP-PEN5);
                    dist(6) = norm(PEP-PEN6);
                    dist(7) = norm(PEP-PEN7);
                    dist(8) = norm(PEP-PEN8);

                    dmax = max(dist);

                    if dmax < Pd          %if width is narrow enough

                        % --------------- Contrast pruning ----------------
                        %can be separated into a function with
                        %(EN1,EN2,EN3,EN4,EN5,EN6,EN7,EN8,P,Pfac,workingimage)
                        %as arguments
                        %output would be (and is) the maximum contrast
                        %level between the neighbors and the pixel
                        %pseudocode
                        %c = max(intensity(q_n)/intensity(p))

                        theta1 = atan2(EN1(2),EN1(1));
                        theta2 = atan2(EN2(2),EN2(1));
                        theta3 = atan2(EN3(2),EN3(1));
                        theta4 = atan2(EN4(2),EN4(1));
                        theta5 = atan2(EN5(2),EN5(1));
                        theta6 = atan2(EN6(2),EN6(1));
                        theta7 = atan2(EN7(2),EN7(1));
                        theta8 = atan2(EN8(2),EN8(1));

                        QN1(1) = round(EN1(1) + Pfac*cos(theta1));
                        QN1(2) = round(EN1(2) + Pfac*sin(theta1));
                        QN2(1) = round(EN2(1) + Pfac*cos(theta2));
                        QN2(2) = round(EN2(2) + Pfac*sin(theta2));
                        QN3(1) = round(EN3(1) + Pfac*cos(theta3));
                        QN3(2) = round(EN3(2) + Pfac*sin(theta3));
                        QN4(1) = round(EN4(1) + Pfac*cos(theta4));
                        QN4(2) = round(EN4(2) + Pfac*sin(theta4));
                        QN5(1) = round(EN5(1) + Pfac*cos(theta5));
                        QN5(2) = round(EN5(2) + Pfac*sin(theta5));
                        QN6(1) = round(EN6(1) + Pfac*cos(theta6));
                        QN6(2) = round(EN6(2) + Pfac*sin(theta6));
                        QN7(1) = round(EN7(1) + Pfac*cos(theta7));
                        QN7(2) = round(EN7(2) + Pfac*sin(theta7));
                        QN8(1) = round(EN8(1) + Pfac*cos(theta8));
                        QN8(2) = round(EN8(2) + Pfac*sin(theta8));

                        %finds the luminance value for the points
                        Pc = greenimage(sub2ind(size(greenimage),P(1),P(2)));
                        QN1c = greenimage(sub2ind(size(greenimage),QN1(1),QN1(2)));
                        QN2c = greenimage(sub2ind(size(greenimage),QN2(1),QN2(2)));
                        QN3c = greenimage(sub2ind(size(greenimage),QN3(1),QN3(2)));
                        QN4c = greenimage(sub2ind(size(greenimage),QN4(1),QN4(2)));
                        QN5c = greenimage(sub2ind(size(greenimage),QN5(1),QN5(2)));
                        QN6c = greenimage(sub2ind(size(greenimage),QN6(1),QN6(2)));
                        QN7c = greenimage(sub2ind(size(greenimage),QN7(1),QN7(2)));
                        QN8c = greenimage(sub2ind(size(greenimage),QN8(1),QN8(2)));

                        Qc(1) = QN1c;
                        Qc(2) = QN2c;
                        Qc(3) = QN3c;
                        Qc(4) = QN4c;
                        Qc(5) = QN5c;
                        Qc(6) = QN6c;
                        Qc(7) = QN7c;
                        Qc(8) = QN8c;
                        %Type conversion to get a decimal value, otherwise
                        %C ends up as an integer
                        Qc = double(Qc);
                        C = Qc./double(Pc);
                        Cmax = max(C);
                        
                        if Cmax > Pcontrast     %if the contrast level passes the threshold for this
                            %if all these tests are passed, contrastimage
                            %is lit up with this pixel.
                            %this was made into a 3D matrix as cumulative
                            %errors were found if this wasn't separated

                            contrastimage(rowdex,coldex,layer) = 255;
                            
                        end %end contrast candidtates
                    end     %end distance (width) candidates
                end         %end angle candidates
            end             %end threshold candidates
        end                 %end cols
    end                     %end rows
    
    % -------------------------Length Pruning------------------------------
    %this section finds all the pixels that are lit up in contrastimage and
    %finds the length of the connected chains. If the chain lengths are below
    %Psize, they are removed from the final image.
    
    %The below assigning doesn't work on multi-page matricies, so we
    %separate the layer
    lengthimage = contrastimage(:,:,layer);
    
    %finds all connected pixels with N8 connectivity
    CC = bwconncomp(lengthimage,8);
    
    %finds length of all connected pixels    
    numlength = cellfun(@length,CC.PixelIdxList);
    
    %finds which connected items are larger than Psize
    numlarger = numlength > Psize;
    
    %finds linear index of all connected items less than length Psize
    lessthanindex = find(not(numlarger));
    
    %eliminates all with length less than Pzise
    for dex=1:length(lessthanindex)
            lengthimage(CC.PixelIdxList{lessthanindex(dex)}) = 0;
    end
    %ORs all the layers together
    SEGMENTED = uint8(lengthimage)|SEGMENTED;
         
end

disp('Sobel Cleaning');
disp(' ');

% end %function end
CLEANED = SEGMENTED;
while 1
        
    %------------------------SOBEL CLEANING------------------------------
    %this section cleans the above segmentation 
    %this section can be called as a function like so:
    %CLEANED = sobel(segmented,greenimage);
    %Line 300 will need to be commented as well. 
    %these Sobel kernels and operators could be placed outside of this
    %while loop, but all our results were made with it inside (as we had
    %each part as a separate function at that point) so we're going to
    %leave it as is JUST IN CASE
    
    %sobel kernels
    G1 = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
    G2 = [1, 2, 1; 0, 0, 0; -1, -2, -1];

    %sobel operators
    Gx = conv2(double(greenimage),double(G1),'same');
    Gy = conv2(double(greenimage),double(G2),'same');
    G = hypot(Gx,Gy);
    Gorient = atan2d(Gy,Gx);
    %end of "stuff that could be outside the loop"
    
    for mindex = 2:rowmax-1
        for nindex = 2:colmax-1

            %these are the neighbour locations, orientations, and gradients of
            %each of the neighbouring pixels 

            %top left
            neighbours(1) = CLEANED(mindex-1, nindex-1);
            orient(1) = Gorient(mindex-1, nindex-1);
            grad(1) = G(mindex-1, nindex-1);

            %top
            neighbours(2) = CLEANED(mindex-1, nindex);
            orient(2) = Gorient(mindex-1, nindex);
            grad(2) = G(mindex-1, nindex);

            %top right
            neighbours(3) = CLEANED(mindex-1, nindex+1);
            orient(3) = Gorient(mindex-1, nindex+1);
            grad(3) = G(mindex-1, nindex+1);

            %left
            neighbours(4) = CLEANED(mindex, nindex-1);
            orient(4) = Gorient(mindex, nindex-1);
            grad(4) = G(mindex, nindex-1);

            %right
            neighbours(5) = CLEANED(mindex, nindex+1);
            orient(5) = Gorient(mindex, nindex+1);
            grad(5) = G(mindex, nindex+1);

            %bottom left
            neighbours(6) = CLEANED(mindex+1, nindex-1);
            orient(6) = Gorient(mindex+1, nindex-1);
            grad(6) = G(mindex+1, nindex-1);

            %bottom
            neighbours(7) = CLEANED(mindex+1, nindex);
            orient(7) = Gorient(mindex+1, nindex);
            grad(7) = G(mindex+1, nindex);

            %bottom right
            neighbours(8) = CLEANED(mindex+1, nindex+1);
            orient(8) = Gorient(mindex+1, nindex+1);
            grad(8) = G(mindex+1, nindex+1);

            if CLEANED(mindex, nindex) > 0    %if vessel candidate
                for count = 1:8
                    if neighbours(count) == 0       %if adjacent to background
                        if abs(orient(count) - Gorient(mindex,nindex)) < .1  %if orientation of grad is same
                           if grad(count) > G(mindex,nindex)     %if 
                                CLEANED(mindex, nindex) = 0;
                            end
                        end                    
                    end
                end
            end
        end
    end
    %tests to ensure each iteration actually thins the vessels. 
    %Returns a Boolean value    
    test = xor(CLEANED,SEGMENTED);
    %checks if any changes are made, concatenates the above tester to a
    %single value
    maxi = max(max(test));
    %if no change is detected
    if maxi == 0
        break;
    else
        SEGMENTED = CLEANED;
    end
%end cleaning section
end %while end

%time taken to segment image
ender = toc(starter);
mins = floor(ender/60);
secs = mod(ender,60);

disp('Time taken for image segmentation')
ddisp = sprintf('%d:%d minutes', mins,secs);
disp(ddisp);
disp('Thank you for your patience');
%displays final segmented image
imtool(CLEANED);
%creates filename for segmented image to write to file
OUTNAME = strcat(INNAME,'_segmented.png');
%writes CLEANED to file
imwrite(CLEANED,OUTNAME);

end %function end