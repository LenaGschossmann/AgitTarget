function spacePress(spaceKey)
    enterStroke = 0 ;
    KbWait([],1);
%     keyCode(:)= 0;
    while(enterStroke == 0)
%         [secs, Time, keyCode] = KbWait;
        [kDown, secs, keyCode, ds] = KbCheck;
        if (find(keyCode == 1) == spaceKey)
            enterStroke = 1;
        end
    end
    enterStroke = 0;
end