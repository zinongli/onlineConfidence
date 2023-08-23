% clear all
Screen('Preference', 'SkipSyncTests', 1); 
cd('C:\Users\labadmin\Documents\onlineConfExperiment');

subj = 'pilot';  
dateTime = clock;                %get  s time for seed  
rng(sum(100*dateTime) );
expName = 'practice';
session = 01;
redoCalib = 0;

[displayInfo] = startExp(subj,datetime,rng);
[displayInfo] = screenVisuals(displayInfo);
if exist(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat']) && redoCalib == 0
    load(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat']) %load calibration incase of restart
    load(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_calibration.mat'])
    
else
    [tform, calibration,startPhase] = penCalib(displayInfo);
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_tform.mat'],'tform')
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_calibration.mat'],'calibration')
end
mode = 0; % lift = 1, slide = 0
%%
start_size = 20;
cursor_size = 5;
pixellength = 0.248;
wait = 1;
patience = 0.5;
topBuff = [0 0 displayInfo.screenXpixels displayInfo.screenAdj/2]; %black bar at top of screen
bottomBuff = [0 displayInfo.screenYpixels-displayInfo.screenAdj/2 displayInfo.screenXpixels displayInfo.screenYpixels]; %black bar at bottom of screen

%% Task Parameters
dists_n = 3;
proj2tablet = 1.886;
pixellength = 0.248;
proj2mm = proj2tablet .* pixellength;
UniRandRadius = 70 .* proj2mm;
edgesize = 50;
sizes_n = 6;
WindowWidth = 1024; % projector window width
yCenter = 384; % projector screen y center
rep = 10;
totalTime = 3;
block_n = 6;

distances = linspace(edgesize,WindowWidth-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,sizes_n*rep) .* proj2mm;

% size granularity TBD, 5:10:55 for now
target_sizes = (5:10:55) ./ pixellength ./ proj2tablet;
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';

lifespan = repmat(totalTime,1,block_n);
%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
traXtotal = [];
traYtotal = [];
testtimes = zeros(1,10000); % 10 seconds
framerate = Screen('NominalFrameRate',displayInfo.window);
frames = framerate * 5; % start/preparing page time out at 5 seconds
instruct = 'Good luck, try hard, and have fun!';
HideCursor;
Screen('FillRect', displayInfo.window, displayInfo.blackVal);
 
while true
    DrawFormattedText(displayInfo.window,instruct,'center','center', displayInfo.whiteVal); 
    Screen('Flip', displayInfo.window);
    [~,~,b] = GetMouse;
    if b(1)
        break
    end
    if KbCheck
        Screen('CloseAll')
        break
    end
end

for j = 1:block_n
    seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
    randdists = distances(seeds);
    randdists = randdists(:);
    randsizes = target_sizes(seeds);
    randsizes = randsizes(:);
    params = NaN(length(randdists),11);
    trax = NaN(length(randdists),round(framerate * (wait+lifespan(j)+patience)));
    tray = NaN(length(randdists),round(framerate * (wait+lifespan(j)+patience)));
    trial_n = length(randdists);
    trials = ones(1,trial_n);
    i = 0;
    DrawFormattedText(displayInfo.window,['Next Block: ' num2str(lifespan(j)) ' seconds interval'],'center','center',displayInfo.whiteVal); % not sure how to get this centered yet
    Screen('Flip', displayInfo.window);
    pause(2);
    while sum(trials) > 0
        i = i+1;
        stage = 0;
        frame = 0;
        if i == trial_n + 1
            randdists = [randdists ; randdists(trials==true,:)];
            randsizes = [randsizes ; randsizes(trials==true,:)];
            params = [params ; params(trials==true,:)];
            trax = [trax ; trax(trials==true,:)];
            tray = [tray ; tray(trials==true,:)];
            wrong_n = sum(trials);
            origin_trial_n = trial_n;
            trial_n = trial_n + wrong_n;
            trials = zeros(1,trial_n);
            trials(1,origin_trial_n+1:end) = 1;
        end
        if stage == 0
            while true
                [~,~,keyCode] = KbCheck;
                if find(keyCode) == 27 % KbName(27) = 'ESCAPE'
                    Screen('CloseAll');
                    ShowCursor;
                    break
                end
                [x,y,buttons] = GetMouse(displayInfo.window2);
                if rem(i,2)
                    startpos = [displayInfo.windowRect(3)-edgesize,displayInfo.yCenter];
                    theta = -pi/12 + (pi/6) * rand(1);
                    rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);
                    params(i,1:2) = startpos - offset;
                else
                    startpos =  [edgesize,displayInfo.yCenter];
                    theta = -pi/12 + (pi/6) * rand(1);
                    rho = randdists(i)-UniRandRadius + UniRandRadius * 2 * rand(1);
                    [offset(1),offset(2)] = pol2cart(theta,rho);
                    params(i,1:2) = startpos + offset;
                end
%                 params(i,10) = randsizes(i);
%                 switch_size = switch_scale * params(i,10);
                Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 1],[],1);
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                if buttons(1)
                    if sqrt(sum((xy - startpos).^2)) <= start_size % if in start area
                        stage = 1;
                        break
                    end
                end
                Screen('Flip', displayInfo.window);
            end
        end
        % test to get the effective refresh rate ( the rate at which the
        % Screen() flips )
        %     testtimes = nonzeros(testtimes);
        %     ts = [0;testtimes];
        %     t_diff = testtimes - ts(1:end-1);
        %     plot(t_diff)
        %     rfrate = 1/mean(t_diff);
        %     worst_rfrate = 1/max(t_diff);
        if stage == 1
            Screen('FillRect', displayInfo.window, displayInfo.blackVal);
            Screen('Flip', displayInfo.window);
            time = GetSecs;
            story = [wait, wait+lifespan(j), wait+lifespan(j)+patience];
            t = 1;
            onset_recorded = 0;
            switch_recorded = 0;
            [x,y,~] = GetMouse(displayInfo.window2);
            tSize = params(i,10);
            for frame = 1: framerate * (story(3)+2) 
                cache = [x,y];
                [x, y, buttons] = GetMouse(displayInfo.window2);
                locdiff = sqrt(sum((cache - [x,y]).^2));
                trax(i,frame) = x;
                tray(i,frame) = y;
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                
                if frame <= framerate * (story(3)+2)
                    Screen('DrawDots', displayInfo.window, startpos, start_size, [1 1 1],[],1);
                    if frame >= framerate * story(1)
                        percent_score = max(1-((frame ./ framerate)-story(1)) / lifespan(j),0);
                        Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,[0 1 0],[],1);
%                         Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
                        Screen('DrawLine', displayInfo.window, [1 1 1], displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                    end
                    Screen('Flip', displayInfo.window);
                    if frame > framerate * story(3) 
                        DrawFormattedText(displayInfo.window,'Too Slow!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                        Screen('Flip',displayInfo.window);
                        trials(i) = 1;
                        params(i,4) = NaN; 
                        pause(1)
                        break
                    end
                    if norm(xy - startpos) >= start_size
                        if buttons(1)+mode ==0
                            DrawFormattedText(displayInfo.window,'Stylus Lifted!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                            Screen('Flip', displayInfo.window);
                            trials(i) = 1;
                            pause(1)
                            break   
                        end
                        if frame < framerate * story(1)
                            DrawFormattedText(displayInfo.window,'Too Early!','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
                            Screen('Flip', displayInfo.window);
                            trials(i) = 1;
                            pause(1)
                            break
                        elseif ~onset_recorded
                            onset_t = frame / framerate;
%                             abs_onset_time = GetSecs;
                            params(i,4) = onset_t;
                            onset_recorded = 1;
                        end
%                         [keyIsDown,~,keyCode] = KbCheck;  % for detecting
%                         keystroke switch
                        
                        if (locdiff <= speedthreshold/framerate && ~mode) || (buttons(1) && mode)
                            if norm(xy - startpos) < randdists(i)/2
                                DrawFormattedText(displayInfo.window,'Not Even Close :(','center','center', displayInfo.whiteVal);
                                Screen('Flip', displayInfo.window);
                                trials(i) = 1;
                                pause(1)  
                            else
                                end_t = frame / framerate;
                                endpos = [x y];
                                hit = 1;
                                bar_color = [1 1 1];
                                if norm(xy - params(i,1:2)) > tSize/2
                                    hit = 0;
                                    bar_color = [1 0 0];
                                end
                                if switch_recorded && hit
                                    end_t = end_t + penalty;
                                    bar_color = [1 1 0];
                                    for k = 1:framerate * penalty
                                        percent_score = max(1-(((frame+k) ./ framerate)-story(1)) / lifespan(j),0);
                                        Screen('DrawDots',displayInfo.window, params(i,1:2), switch_size,[0 1 0],[],1);
                                        Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                        Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                        DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)) '/' num2str(length(distances)) ' finished'],'center','center', displayInfo.whiteVal);
                                        Screen('Flip', displayInfo.window);
                                    end
                                    score = percent_score * 10 * hit;
                                else
                                    Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,[0 1 0],[],1);
%                                     Screen('FrameOval',displayInfo.window, [1 1 0], [params(i,1)-switch_size./2,params(i,2)-switch_size./2,params(i,1)+switch_size./2,params(i,2)+switch_size./2]);
                                    Screen('DrawDots', displayInfo.window, xy, 5, [1 0 0],[],1);
                                    Screen('DrawLine', displayInfo.window, bar_color, displayInfo.xCenter - percent_score * scorebar_length,displayInfo.yCenter-200, displayInfo.xCenter + percent_score * scorebar_length,displayInfo.yCenter-200,5);
                                    if hit
                                        DrawFormattedText(displayInfo.window,['Score = ' num2str(percent_score * 10)],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    else
                                        DrawFormattedText(displayInfo.window,['Miss :('],'center',displayInfo.yCenter-220, displayInfo.whiteVal);
                                    end
                                    DrawFormattedText(displayInfo.window,[num2str(length(trials)-sum(trials)+1) '/' num2str(length(distances)*2) ' finished'],'center','center', displayInfo.whiteVal);
                                    Screen('Flip', displayInfo.window);
                                    score = percent_score * 10 * hit;
                                end
                                rest_of_trial = story(3) - end_t;
                                pause(rest_of_trial);
                                params(i,3) = 0;
                                params(i,5) = end_t;
                                params(i,6:7) = endpos;
                                params(i,8:9) = startpos;
                                params(i,10) = randsizes(i);
                                params(i,11) = score;
                                trials(i) = 0;
                            end
                            break
                        end                            
                    end
                end
%                 testtimes(t) = GetSecs; % for testing temporal resolution
%                 t = t+1;
                [~,~,keyCode] = KbCheck;
                if find(keyCode) == 27 % KbName(27) = 'ESCAPE'
                    Screen('CloseAll');
                    ShowCursor;
                    break
                end
            end
        end
    end
    data = [data;params];
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_rawtotal.mat'],'data');
    xTrajTelomere = NaN(size(trax,1),size(traXtotal,2));
    xTrajTelomere(1:size(trax,1),1:size(trax,2)) = trax;
    traXtotal = [traXtotal;xTrajTelomere];
    yTrajTelomere = NaN(size(tray,1),size(traYtotal,2));
    yTrajTelomere(1:size(tray,1),1:size(tray,2)) = tray;
    traYtotal = [traYtotal;yTrajTelomere];
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_traXtotal.mat'],'traXtotal')
    save(['data_onlineConf\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_traYtotal.mat'],'traYtotal')
    while true
        DrawFormattedText(displayInfo.window,'Block finished. Press any key to proceed to next block.','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
        Screen('Flip', displayInfo.window);
        if KbCheck
            break
        end
    end
end
Screen('CloseAll');
ShowCursor;