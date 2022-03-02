function Fastsavetiff(filename,images)

    buffer = uint16(images);
    fTIF = Fast_BigTiff_Write(filename,1,0);
%     msg = 0;
    for ct = 1:size(buffer,3)
%         fprintf(1,repmat('\b',[1,msg]));
%         msg = fprintf(1,'%.0f/%.0f',ct,75000);
        %noise = uint16(2000+randn(size(IM))*1000);
        %fTIF.WriteIMG(IM+noise);
        fTIF.WriteIMG(buffer(:,:,ct)');
    end
    fTIF.close;
    disp('save tiff finished')
    % fprintf(1,repmat('\b',[1,msg]));
    % fprintf(1,'\nWrite %.0f bytes in %.0f seconds\n',B*N,t);
    % fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));