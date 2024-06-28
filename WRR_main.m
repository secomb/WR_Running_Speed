%
%   Code to accompany manuscript "An empirical model for world record running speeds with distance, age, and sex: anaerobic and aerobic contributions to performance"
%   by Tuhin K. Roy [roy.tk(at)mayo.edu], Michael J. Joyner, Jonathon W. Senefeld, Chad C. Wiggins, Timothy W. Secomb (under review)
%   May 2024
%   Dependencies: rmsecalc.m, sigage5.m
%   Data file: Running4.xlsx


close all
clear all
clc
format compact

RM = xlsread('Running4.xlsx','Men Outdoor - All Info');
RF = xlsread('Running4.xlsx','Women Outdoor - All Info');

distcol = 1;
timecol = 3;
dm = RM(:,distcol);
tm = RM(:,timecol);
spm = dm./tm;                               % Speed, males
df = RF(:,distcol);
tf = RF(:,timecol);
spf = df./tf;                               % Speed, females

% Split data from athletefirst.com
splitm = [5.47 9.58; 9.92 19.19];           % 100m and 200m, Bolt 2009 IAAF WC
splitf = [5.97 10.49; 11.18 21.34];         % 100m Griffith Joyner 1988 USA Olympic Trials; 200m Griffith-Joyner 1988 Olympic Games
splitspm = [50/(splitm(1,2)-splitm(1,1)) 100/(splitm(2,2)-splitm(2,1))];
splitspf = [50/(splitf(1,2)-splitf(1,1)) 100/(splitf(2,2)-splitf(2,1))];

% Set up Figure 10 for tiled layout
f10 = figure(10);
f10.Position = [1, 1, 1000, 600];
t10 = tiledlayout(1,2);

% Set up Figure 15 for tiled layout
f15 = figure(15);
f15.Position = [1, 1, 1400, 600];
t15 = tiledlayout(1,2);

% Figures: Men/red on left, Women/blue on right
for gender = 1:2                                            % 1 = Men, 2 = Women
    if (gender == 1)
        distlen = 11;
        distcat = RM(1:distlen);
        agecat = 30:5:105;                                  % Lower limits of age categories [<35 35-40 ... 105-110]
        agelen = length(agecat);
        tmat = reshape(tm,[distlen,agelen]);
        dmat = reshape(dm,[distlen,agelen]);
        spmatm = reshape(spm,[distlen,agelen]);
    else
        distcat = RF(1:distlen);
        agecat = 30:5:105;
        agelen = length(agecat);
        tmat = reshape(tf,[distlen,agelen]);
        dmat = reshape(df,[distlen,agelen]);
        spmat = reshape(spf,[distlen,agelen]);
        agecat = 30:5:100;                                  % Lower limits of age categories [<35 35-40 ... 100-105]
        agelen = length(agecat);
        tmat = tmat(1:end,1:agelen);
        dmat = dmat(1:end,1:agelen);
        spmatf = spmat(1:end,1:agelen);
    end
    
    if (gender == 1)
        spmat = spmatm;
    else
        spmat = spmatf;
    end
    
    spmatcorr = spmat;
    agefac = 1;
    age = 1;
    if (gender == 1)
        spmatcorr(1,age) = spmat(1,age)*(1+((splitspm(1)/spmat(1,1))-1)*agefac);
        spmatcorr(2,age) = spmat(2,age)*(1+((splitspm(2)/spmat(2,1))-1)*agefac);
    else
        spmatcorr(1,age) = spmat(1,age)*(1+((splitspf(1)/spmat(1,1))-1)*agefac);
        spmatcorr(2,age) = spmat(2,age)*(1+((splitspf(2)/spmat(2,1))-1)*agefac);
    end
    
    % Distance fit for WR (youngest) group
    
    xval = dmat(:,1);
    agemax = agelen;
    colordata = jet(agemax);
    coeffs = [];
    
    ageindex = 1;
    normyval = spmatcorr(:,ageindex)/spmatcorr(1,ageindex);
    yval = spmatcorr(:,ageindex);
    % semilogx(xval,yval,'ro')
    
    figure(3)
    dref = 100;
    eqdist1 = @(Vref,B,b,S,d,x)(Vref - S*log10(x/dref) + B*(1+(dref/d).^b)./(1+(x/d).^b));
    distfitoptions = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1,1,1,1,100]);
    [distfit1, distgof1, distgofout] = fit(xval,yval,eqdist1,distfitoptions)
    [[feval(distfit1,xval)] yval]';
    [distrmserr, distrmsalt] = rmsecalc([feval(distfit1,xval)],yval,length(coeffvalues(distfit1)))
    % [distfit1, distgof1] = fit(xval,yval,eqdist1,'StartPoint',[1,1,1,1,100])
    % coeffs = [coeffs' coeffvalues(distfit1)']';
    Vref0 = distfit1.Vref
    B0 = distfit1.B
    b0 = distfit1.b
    S0 = distfit1.S
    d0 = distfit1.d
    if (gender == 1)
        Vref0_distfit_men = Vref0;
        B0_distfit_men = B0;
        b0_distfit_men = b0;
        S0_distfit_men = S0;
        d0_distfit_men = d0;
        distrmserr_men = distrmserr;
    else
        Vref0_distfit_women = Vref0;
        B0_distfit_women = B0;
        b0_distfit_women = b0;
        S0_distfit_women = S0;
        d0_distfit_women = d0;
        distrmserr_women = distrmserr;
    end
    xrange = 10.^[log10(min(xval)):0.01:log10(max(xval))];
    subplot(1,2,gender)
    if (gender == 1)
        fitline = semilogx(xrange,feval(distfit1,xrange),'LineStyle','-','Color','b');
    else
        fitline = semilogx(xrange,feval(distfit1,xrange),'LineStyle','-','Color','r');
    end
    % fitline = semilogx(xrange,feval(distfit1,xrange),'Color',colordata(ageindex,:),'LineStyle','-');
    % color = get(fitline, 'Color');
    grid on
    hold all
    % semilogx(xval,yval,'s'),'Color',color)
    if (gender == 1)
        semilogx(xval,yval,'sb')
    else
        semilogx(xval,yval,'sr')
    end
    hold on
    yaer = Vref0 - S0.*log10(xval/dref);
    if (gender == 1)    
        semilogx(xval,yaer,'b')
    else
        semilogx(xval,yaer,'r')
    end
    yvalorig = spmat(:,ageindex);
    hold on
    plot(xval(1:2),yvalorig(1:2),'^')   % Actual uncorrected value
    plot(xval(1:3),yvalorig(1:3),'-')   % Actual uncorrected value
    ylabel('Speed (m/s)')
    xlabel('Distance (m)')
    ylim([5 13])
    if (gender == 1)
        title('A')
        subplot(1,2,2)
        fitline = semilogx(xrange,feval(distfit1,xrange),'Color',colordata(ageindex,:),'LineStyle','--');
        semilogx(xval,yaer,'b')
        hold on
    else       
        semilogx(xval,yaer,'r')
        title('B')
    end
    
    figure(9)
    subplot(1,2,1)
    xlabel('Distance (m)')
    boostspeed = @(x)(distfit1.B*(1+(dref/distfit1.d).^distfit1.b)./(1+(x/distfit1.d).^distfit1.b));
    if (gender == 1)
        semilogx(xrange,feval(boostspeed,xrange),'Color','b')
    else
        semilogx(xrange,feval(boostspeed,xrange),'Color','r')
    end
    % axis([-Inf Inf 0 boostspeed(min(xval))])
    ylabel('Increase in running speed (m/s)')
    title('A')
    hold on;
    
   % xinp = input("end of plot");
    subplot(1,2,2)
    timeadv = @(x)(x./(distfit1.Vref - distfit1.S.*log10(x./dref)) - x./feval(distfit1,x)');
    if (gender == 1)
        semilogx(xrange,feval(timeadv,xrange),'Color','b')
    else
        semilogx(xrange,feval(timeadv,xrange),'Color','r')
    end
    %    semilogx(xrange,xrange./(agecoeff(agecat(1))*(fitcurve5.Vref - fitcurve5.S.*log10(xrange/dref)))-xrange./(agecoeff(agecat(1))*feval(eqdist5,xrange)))
    ylabel('Decrease in running time (s)')
    xlabel('Distance (m)')
    title('B')
    hold on;
    % xinp = input("end of plot");
    
    iter = 0;
    done = false;
    while (~done)
        spmatcorr = spmat;
        iter = iter + 1
        maxage = 105;
        for age = 1:agelen
            agefac = (maxage-agecat(age))/(maxage-agecat(1));
            % Correct speed for startup time by using ratio of split speed to avg for 100m and 200m
            % Assume linear decrease in effect of startup time, with factor
            % ranging from 1 for youngest group to 0 for age 105
            if (gender == 1)
                spmatcorr(1,age) = spmat(1,age)*(1+((splitspm(1)/spmat(1,1))-1)*agefac);
                spmatcorr(2,age) = spmat(2,age)*(1+((splitspm(2)/spmat(2,1))-1)*agefac);
            else
                spmatcorr(1,age) = spmat(1,age)*(1+((splitspf(1)/spmat(1,1))-1)*agefac);
                spmatcorr(2,age) = spmat(2,age)*(1+((splitspf(2)/spmat(2,1))-1)*agefac);
            end
        end
        
        figure(2);
        xval = agecat;
        yval = spmatcorr(1,:);
        
        x0 = xval(1);
        y0 = yval(1);
        
        % Age fit of corrected speeds for WR data to find A0 and Yref0
        
        ageeqn = @(A,C,Yref,x) (C*(1-A*exp(x/Yref)))
        agefitoptions = fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[1,y0,x0])
        [agefit, agegof, ageout] = fit(xval',yval',ageeqn,agefitoptions)
        [agermserr, agermsalt] = rmsecalc(feval(agefit,xval),yval',length(coeffvalues(agefit)))
        
        A0 = agefit.A;
        C0 = agefit.C;
        Yref0 = agefit.Yref;
        
        subplot(1,2,gender)
        if (gender == 1)
            subplot(1,2,gender)
            plot(agefit,'b')
            hold on
            plot(xval,yval,'ob','MarkerSize',4)
            hold off
            agefitm = agefit;
            A0_agefit_men = A0;
            C0_agefit_men = C0;
            Yref0_agefit_men = Yref0;
            agermserr_men = agermserr;
            legend off
            title('A')
            xlim([20 105])
            ylim([0 13])
        else
            plot(agefit,'r')
            hold on
            plot(xval,yval,'or','MarkerSize',4)
            hold off
            agefitf = agefit;
            A0_agefit_women = A0;
            C0_agefit_women = C0;
            Yref0_agefit_women = Yref0;
            agermserr_women = agermserr;
            legend off
            title('B')
            xlim([20 105])
            ylim([0 13])
        end
        xlabel('Age (y)')
        ylabel('Speed (m/s)')
        % [feval(agefit,xval) yval']'
        
        k0 = 0;
        %        figure(8)
        xval = dmat(1:distlen)';
        yval = agecat';
        zval = spmatcorr;
        [X,Y] = meshgrid(xval,yval);
        %        surf(log(X),Y,zval')
        [x0,y0,z0] = prepareSurfaceData(xval,yval,zval');
        
        warnmsg = warning('query','last');
        warnid = warnmsg.identifier;
        warning('off',warnid);
        
        % Fit assuming that anaerobic boost is independent of age
        
        % Use range of values for youngest group category
        
        figure(6)
        fitformula = 'sigage5(Xdist, Yage, Vref, B, b, S, d, Yref, A)';
        ft5 = fittype(fitformula,'numindep',2,'independent',{'Xdist','Yage'},'dependent',{'Z'}, ...
            'coefficient',{'Vref','B','b','S','d','Yref','A'});
        
        for youngestage = 18:30
            xval = dmat(1:distlen)';
            yval = agecat';
            yval(1) = youngestage;
            zval = spmatcorr;
            [x0,y0,z0] = prepareSurfaceData(xval,yval,zval');
            ft5options = fitoptions('Method','NonlinearLeastSquares','StartPoint',[Vref0,B0,b0,S0,d0,Yref0,A0]);
            [fitcurve5, gof5] = fit([x0,y0],z0,ft5,ft5options);
            
            gofage(youngestage) = gof5.rmse;
            
            eqdist5 = @(x)(fitcurve5.Vref - fitcurve5.S*log10(x/dref) + fitcurve5.B*(1+(dref/fitcurve5.d).^fitcurve5.b)./(1+(x/fitcurve5.d).^fitcurve5.b));
            agecoeff = @(y)(1-fitcurve5.A*exp(y/fitcurve5.Yref));
            
            youngestage;
            agecoeff(youngestage);
            ageindex = youngestage - 18 + 1;
            subplot(1,2,gender)
            fitline = semilogx(xrange,agecoeff(youngestage)*feval(eqdist5,xrange),'Color',colordata(ageindex,:),'LineStyle','--');
            color = get(fitline, 'Color');
            hold on
            ageindex = 1;
            yval = spmatcorr(:,ageindex);
            % semilogx(xval,yval,'o','Color',color)
            % Calculate rmse for youngest age group youngestage
            yfit = agecoeff(youngestage)*feval(eqdist5,xval);
            %             fiterrsq = 0;
            %             for i = 1:distlen
            %                 fiterrsq = fiterrsq + (yfit(i) - yval(i))^2;
            %             end
            %             fitrmse(youngestage) = sqrt(fiterrsq/distlen)
            fitrmse(youngestage) = rmsecalc(yfit,yval,0);
            % fitval = semilogx(xval,yfit,'*')
            grid on
        end
        
        fitrmse';
        [bestage,bestageindex] = min(fitrmse(18:30));
        bestage
        youngestageref = bestageindex + 18 - 1
        if (youngestageref == agecat(1))
            done = true;
            figure(7)
            subplot(1,2,2)
            if gender == 1
                plot(18:30,fitrmse(18:30),'b')
                hold on
                plot(bestageindex + 18 - 1,bestage,'sb')
            else
                plot(18:30,fitrmse(18:30),'r')
                plot(bestageindex + 18 - 1,bestage,'sr')
                hold off
            end
            xlabel('Assigned WR Cohort Age')
            ylabel('RMS deviation (m/s), WR group')
            title('B')
            subplot(1,2,1)
            if gender == 1
                plot(18:30,gofage(18:30),'b')
                hold on
                plot(bestageindex + 18 - 1,gofage(bestageindex + 18 -1),'sb')
            else
                plot(18:30,gofage(18:30),'r')
                plot(bestageindex + 18 - 1,gofage(bestageindex + 18 -1),'sr')
                hold off
            end
            xlabel('Assigned WR Cohort Age')
            ylabel('RMS deviation (m/s), all ages')
            title('A')
        end
        agecat(1) = youngestageref
    end
    %   figure
    xval = dmat(1:distlen)';
    yval = agecat';
    zval = spmatcorr;
    [x0,y0,z0] = prepareSurfaceData(xval,yval,zval');
    
    % Fit over entire dataset
    
    ft5options = fitoptions('Method','NonlinearLeastSquares','StartPoint',[Vref0,B0,b0,S0,d0,Yref0,A0]);
    [fitcurve5, gof5, fitoutput5] = fit([x0,y0],z0,ft5,ft5options);
    zfit = feval(fitcurve5,[x0,y0]);
    [surfrmserr, surfrmsalt] = rmsecalc(zfit,z0,length(coeffvalues(fitcurve5)));
    Vref0 = fitcurve5.Vref;
    B0 = fitcurve5.B;
    b0 = fitcurve5.b;
    S0 = fitcurve5.S;
    d0 = fitcurve5.d;
    Yref0 = fitcurve5.Yref;
    A0 = fitcurve5.A;
    if (gender == 1)
        Vref0_surffit_men = Vref0;
        B0_surffit_men = B0;
        b0_surffit_men = b0;
        S0_surffit_men = S0;
        d0_surffit_men = d0;
        Yref0_surffit_men = Yref0;
        A0_surffit_men = A0;
        surfrmserr_men = surfrmserr;
    else
        Vref0_surffit_women = Vref0;
        B0_surffit_women = B0;
        b0_surffit_women = b0;
        S0_surffit_women = S0;
        d0_surffit_women = d0;
        Yref0_surffit_women = Yref0;
        A0_surffit_women = A0;
        surfrmserr_women = surfrmserr;
    end
    agemax = agelen;
    colordata = jet(agemax);
    
    ageindex = 1;
    yval = spmatcorr(:,ageindex);
    
    rmsarray = zeros(agelen,1);
    % Make semilog plots by age
    figure(10) % Tiled layout
    lnw = 1;
%    subplot(1,2,gender)
    nexttile(gender);
    for age = 1:agelen
        fitplot = semilogx(xrange,feval(eqdist5,xrange)*agecoeff(agecat(age)),'LineWidth',lnw);
        [rmserr, rmsalt] = rmsecalc(feval(eqdist5,xval)*agecoeff(agecat(age)),zval(:,age),length(coeffvalues(fitcurve5)));
        % [feval(eqdist5,xval)*agecoeff(agecat(age)) zval(:,age)]'
        rmsarray(age) = rmserr;
        color = get(fitplot, 'Color');
        hold on;
        semilogx(xval,zval(:,age),'o','MarkerSize',4,'Color',color);
    end
    agecatstr = string(agecat);
    for i = 1:length(agecatstr)
        agecatleg(i*2-1) = agecatstr(i);
        agecatleg(i*2) = "";
    end
    agecatleg(1) = "WR"
    [leg, legobj, ~, ~] = legend(agecatleg);
    hlegobj = findobj(legobj,'type','line')
    set(hlegobj,'LineWidth',lnw);
    leg.ItemTokenSize = [10,72]; % Default [30,18]
    leg.Box = 'off'
    ylim([0 13])
    xlim([0 10^6])
    if (gender == 1)
        title('A');
    else
        title('B');
    end
    %    title('Predicted speed')
    xlabel('Distance (m)')
    ylabel('Speed (m/s)')
    
    figure(15)
    ax = nexttile(gender)
%    ax = subplot(1,2,gender)
    agecatstr = string(agecat);
    agecatstr(1) = "WR"
    rmsarray'
    if (gender == 1)
        bar(rmsarray',1,'b')
    else
        bar(rmsarray',1,'r')
    end
    set(ax,'XTick',1:length(agecatstr),'XTickLabel',agecatstr,'TickLength',[0.01 0])
    ylim([0 0.7]);
    xlabel('Age (y)')
    if (gender == 1)
        title('A')
        ylabel('RMS deviation (m/s)')
    else
        title('B')
        ylabel('RMS deviation (m/s)')
%         pos = get(ax,'Position');
%         set(ax, 'Units', 'normalized')
%         set(ax, 'Position', pos - [0.05 0 0 0])
%         pos = get(ax,'Position')
    end
end

Table1 = [C0_agefit_men A0_agefit_men Yref0_agefit_men agermserr_men; ...
    C0_agefit_women A0_agefit_women Yref0_agefit_women agermserr_women]'

Table2 = [Vref0_distfit_men S0_distfit_men B0_distfit_men b0_distfit_men d0_distfit_men distrmserr_men; ...
    Vref0_distfit_women S0_distfit_women B0_distfit_women b0_distfit_women d0_distfit_women distrmserr_women]'

Table3 = [Vref0_surffit_men S0_surffit_men B0_surffit_men b0_surffit_men d0_surffit_men A0_surffit_men Yref0_surffit_men surfrmserr_men;
    Vref0_surffit_women S0_surffit_women B0_surffit_women b0_surffit_women d0_surffit_women A0_surffit_women Yref0_surffit_women surfrmserr_women]'

    fig1 = figure(3);
    exportgraphics(fig1,'Fig1.tif','Resolution',600)
    
    fig2 = figure(9);
    exportgraphics(fig2,'Fig2.tif','Resolution',600)
   
    fig3 = figure(2);
    exportgraphics(fig3,'Fig3.tif','Resolution',600)
    
    fig4 = figure(7);
    exportgraphics(fig4,'Fig4.tif','Resolution',600)
   
    fig5 = t10;
    exportgraphics(fig5,'Fig5.tif','Resolution',600)
    
    fig6 = t15;
    exportgraphics(fig6,'Fig6.tif','Resolution',600)