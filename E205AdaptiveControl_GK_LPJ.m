function varargout = E205AdaptiveControl_GK_LPJ(varargin)
% E205ADAPTIVECONTROL_GK_LPJ MATLAB code for E205AdaptiveControl_GK_LPJ.fig
%      E205ADAPTIVECONTROL_GK_LPJ, by itself, creates a new E205ADAPTIVECONTROL_GK_LPJ or raises the existing
%      singleton*.
%
%      H = E205ADAPTIVECONTROL_GK_LPJ returns the handle to a new E205ADAPTIVECONTROL_GK_LPJ or the handle to
%      the existing singleton*.
%
%      E205ADAPTIVECONTROL_GK_LPJ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in E205ADAPTIVECONTROL_GK_LPJ.M with the given input arguments.
%
%      E205ADAPTIVECONTROL_GK_LPJ('Property','Value',...) creates a new E205ADAPTIVECONTROL_GK_LPJ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before E205AdaptiveControl_GK_LPJ_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to E205AdaptiveControl_GK_LPJ_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help E205AdaptiveControl_GK_LPJ

% Last Modified by GUIDE v2.5 16-Dec-2014 19:08:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @E205AdaptiveControl_GK_LPJ_OpeningFcn, ...
                   'gui_OutputFcn',  @E205AdaptiveControl_GK_LPJ_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before E205AdaptiveControl_GK_LPJ is made visible.
function E205AdaptiveControl_GK_LPJ_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to E205AdaptiveControl_GK_LPJ (see VARARGIN)

% Choose default command line output for E205AdaptiveControl_GK_LPJ
handles.output = hObject;

% Update handles structure
handles.Jactual = 2;
handles.Jguess = 1;
handles.initialConditions = [0.5 0.5];
handles.ditherAmplitude = 0.01;
handles.ditherFreq = 100;
handles.timeSpan = 50;
handles.plotSpan = 10;
handles.refSignal = 'Pulse';
handles.refAmplitude = 1;
handles.refFreq = 100;
handles.results = {};

guidata(hObject, handles);

% UIWAIT makes E205AdaptiveControl_GK_LPJ wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% % Plot stability boundary
% axes(handles.Phat);
% hold all
% 
% % Calculate and plot stable boundaries of system
% mu_stable_bound = linspace(0, 0.9, 1000);
% sigma_stable_bound = mu_stable_bound./(1-mu_stable_bound);
% plot(mu_stable_bound, sigma_stable_bound, 'red');
% 
% % Calculate and plot Transition from Focus to Node Type
% mu_focus = linspace(0.8, 0.9, 100);
% sigma_focus = (mu_focus-sqrt(4.*(1-mu_focus)))./(1-mu_focus);
% plot(mu_focus, sigma_focus, 'green');
% 
% % Plot Transcritical Bifurcation Point
% line([1 1],[0 2], 'Color', 'blue');
% 
% %Initialize Parameter Axes
% xlim([0 2]);
% ylim([0 2]);
% title('Parameter Space')
% xlabel('mu (Prey overcrowding and disease factor)');
% ylabel('sigma (First order time delay for predator growth)');
% 
% % Add Region Labels
% text(0.3,1.95,'1') 
% text(0.74,1.95,'2')
% text(0.92,1.95,'3')
% text(1.4,1.95,'4') 
% % legend('Stable','Oscillatory','Predator Threshold)');
% 
% % Initialize Phase Portraits
% axes(handles.Theta);
% hold all
% title('Phase Portrait');
% xlabel ('x_1 (Prey Biomass)')
% ylabel ('x_2 (Predator Biomass)')
% axis auto
% 
% % Initialize Time Plots
% axes(handles.Dhat);
% hold all
% title('Time Evolution');
% xlabel ('Time')
% ylabel ('Predator/Prey Biomass')

% --- Outputs from this function are returned to the command line.
function varargout = E205AdaptiveControl_GK_LPJ_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Updates phase plot when anything changes
function handles = updatePhasePlot(handles)

% 
% % Store current values of mu and sigma into variables
% mu = handles.mu;
% sigma = handles.sigma;
% 
% % Readjust sliders to reflect latest values of mu and sigma
% set(handles.mu_slider,'value', mu);
% set(handles.sigma_slider,'value', sigma);
% 
% % Define and clear axes
% axes(handles.Phat);
% cla(handles.Phat);
% 
% % Add Region Labels
% text(0.3,1.95,'1') 
% text(0.74,1.95,'2')
% text(0.92,1.95,'3')
% text(1.4,1.95,'4') 
% 
% % Calculate and plot stable boundaries of system
% mu_stable_bound = linspace(0, 0.9, 1000);
% sigma_stable_bound = mu_stable_bound./(1-mu_stable_bound);
% plot(mu_stable_bound, sigma_stable_bound, 'red');
% 
% % Calculate and plot Transition from Focus to Node Type
% mu_focus = linspace(0.8, 0.9, 100);
% sigma_focus = (mu_focus-sqrt(4.*(1-mu_focus)))./(1-mu_focus);
% plot(mu_focus, sigma_focus, 'green');
% 
% % Plot Transcritical Bifurcation Point
% line([1 1],[0 2], 'Color', 'blue');
% 
% plot(mu,sigma,'ro')
% 
% % Use presolved fixed point locations
% xstar = [0, 1/mu, 1];
% ystar = [0, 0, 1-mu];
% 
% % Define axes
% axes(handles.Theta)
% 
% % Store old x and y axis limits
% old_xlims = xlim;
% old_ylims = ylim;
% 
% % Clear axes
% cla(handles.Theta);
% 
% % Plot fixed points
% plot(xstar,ystar,'ro')
% 
% 
% 
% % Solve with initial conditions
% tspan=[0 handles.timeSpan];
% x0=handles.initialConditions';
% 
% % Define the system
% dx = @(t,x) [x(1) - x(1)*x(2) -  mu*x(1)^2; x(1)*x(2) - x(2) - sigma*x(2) * (x(1) - x(1)*x(2) -  mu*x(1)^2) ];
% 
% % Solve Timeout problem: http://www.mathworks.com/matlabcentral/newsreader/view_thread/119565
% % 1. define a link to an an event function which will stop calculation
% xoverFcn = @(T, Y) MyEventFunction(T, Y); % defined at bottom of code
% % 2. register this function as an event function
% options = odeset('Events',xoverFcn); 
% % 3. start a stopwatch timer, if you already use one, define a new one: tic(ticID)
% tic;
% [t,x]=ode45(dx,tspan,x0,options);
% 
% % Plot phase results
% plot(x(:,1),x(:,2), 'blue')
% 
% % If "Hold Axis Limits" is checked, reapply old axis limits
% if (get(handles.hold_axis_lims, 'Value')) 
%     axis([old_xlims old_ylims])
% else
%     axis auto
%     axLims = axis;
%     axis(max(0,axLims))
% end
% 
% %Plot time domain results
% axes(handles.Dhat)
% 
% % Store old x and y axis limits
% old_xlims = xlim;
% old_ylims = ylim;
% 
% cla(handles.Dhat);
% hold all
% plot(t,x(:,1), 'green');
% plot(t,x(:,2), 'red');
% legend('Prey','Predator');
% 
% % If "Hold Axis Limits" is checked, reapply old axis limits
% if (get(handles.hold_axis_lims, 'Value')) 
%     axis([old_xlims old_ylims])
% else
%     axis auto
% end
% 
% % Update results table
% % trace was found using jacobian at fixed point [1, 1-mu]
% nodeLoc = [1, 1-handles.mu];
% trace = -mu + sigma*(1-mu);
% det = 1-mu;
% if det < 0 % [1, 1-mu] node is not the stable one, the [1/mu, 0] one is.
%     nodeLoc = [1/handles.mu, 0];
%     str = sprintf('Stable Node at [%.3g, %.3g]', nodeLoc);
%     newResult = {handles.mu, handles.sigma, str, 'N/A', 'N/A'};
% elseif trace > 0 % unstable, limit cycle predicted    
%     w = sqrt(det - (trace/2)^2); %Imaginary part of eigenvalues
%     Tpredicted = 2*pi/w;
%     Tmeasured = measureLimitCycle(t,x(:,1));
%     str = sprintf('Limit Cycle around [%.3g, %.3g]', nodeLoc);
%     newResult = {handles.mu, handles.sigma, str, Tpredicted, Tmeasured};
% else % stable, check whether eigenvalues are real or complex
%     has_real_roots = trace^2 - 4*det;
%     if (has_real_roots > 0)
%         str = sprintf('Stable Node at [%.3g, %.3g]', nodeLoc);
%     else
%         str = sprintf('Stable Focus at [%.3g, %.3g]', nodeLoc);
%     end
%     newResult = {handles.mu, handles.sigma, str, 'N/A', 'N/A'};
% end
% 
% if isempty(handles.results) || ~isequal(handles.results(1,:), newResult)
%     handles.results = [newResult; handles.results];
% end
% set(handles.resultsTable,'Data',handles.results)


% --- Executes on mouse press over axes background.
function Phat_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to stabilityBoundary (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Updates values of a and b (reaction constants)
% coords = get(hObject, 'CurrentPoint');
% newMu = coords(1,1);
% newSigma = coords(1,2);
% 
% % Correct clicking off the screen
% newMu = max(newMu,0);
% newMu = min(newMu,2);
% 
% newSigma = max(newSigma,0);
% newSigma = min(newSigma,2);
% 
% handles.mu = newMu;
% handles.sigma = newSigma;
% mu_str = sprintf('%g', handles.mu);
% sigma_str = sprintf('%g', handles.sigma);
% set(handles.mu_disp, 'String', mu_str);
% set(handles.sigma_disp, 'String', sigma_str);
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)

function mu_slider_Callback(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu_disp as text
%        str2double(get(hObject,'String')) returns contents of mu_disp as a double
% handles.mu = get(hObject,'value');
% mu_str = sprintf('%g', handles.mu);
% set(handles.mu_disp, 'String', mu_str);
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function mu_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_slider_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_disp as text
%        str2double(get(hObject,'String')) returns contents of sigma_disp as a double
% handles.sigma = get(hObject,'value');
% sigma_str = sprintf('%g', handles.sigma);
% set(handles.sigma_disp, 'String', sigma_str);
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function sigma_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mu_disp_Callback(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu_disp as text
%        str2double(get(hObject,'String')) returns contents of mu_disp as a double
newMu = str2double(get(hObject,'String'));
if newMu > 2
    set(hObject, 'String', '2')
    newMu = 2;
elseif newMu < 0
    set(hObject, 'String', '0')
    newMu = 0;
end
handles.mu = newMu;
handles = updatePhasePlot(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function mu_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_disp_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_disp as text
%        str2double(get(hObject,'String')) returns contents of sigma_disp as a double
newSigma = str2double(get(hObject,'String'));
if newSigma > 2
    set(hObject, 'String', '2')
    newSigma = 2;
elseif newSigma < 0
    set(hObject, 'String', '0')
    newSigma = 0;
end
handles.sigma = newSigma;
handles = updatePhasePlot(handles);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function sigma_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initCondDisp_Callback(hObject, eventdata, handles)
% hObject    handle to initCondDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initCondDisp as text
%        str2double(get(hObject,'String')) returns contents of initCondDisp as a double

% --- Executes during object creation, after setting all properties.
function initCondDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initCondDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in initCondEdit.
function initCondEdit_Callback(hObject, eventdata, handles)
% hObject    handle to initCondEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% prompt = {'Intial Prey Biomass','Initial Predator Biomass'};
% dlg_title = 'Initial Conditions';
% num_lines = 1;
% ic = handles.initialConditions;
% def = {num2str(ic(1)), num2str(ic(2))};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% assignin('base', 'answer', answer)
% if isempty(answer)
%     %do nothing
% else
%     % Update initial conditions
%     handles.initialConditions = min(max([eval(answer{1}), eval(answer{2})],0),1);
%     ic = handles.initialConditions;
%     icstr = sprintf('Initial Conditions: [%g %g]', ic(1), ic(2));
%     set(handles.initCondDisp, 'String', icstr);
% end
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)


function timeSpanDisp_Callback(hObject, eventdata, handles)
% hObject    handle to timeSpanDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeSpanDisp as text
%        str2double(get(hObject,'String')) returns contents of timeSpanDisp as a double


% --- Executes during object creation, after setting all properties.
function timeSpanDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSpanDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in timeSpanButton.
function timeSpanButton_Callback(hObject, eventdata, handles)
% hObject    handle to timeSpanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% prompt = {'Time Span (s)'};
% dlg_title = 'Time Span';
% num_lines = 1;
% tspan = handles.timeSpan;
% def = {num2str(tspan)};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% assignin('base', 'answer', answer)
% if isempty(answer)
%     %do nothing
% else
%     % Update time span
%     handles.timeSpan = min(max(eval(answer{1}),1),300);
%     tspan = handles.timeSpan;
%     tstr = sprintf('Time Span: %g', tspan);
%     set(handles.timeSpanDisp, 'String', tstr);
% end
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)

% --- Executes on button press in hold_axis_lims.
function hold_axis_lims_Callback(hObject, eventdata, handles)
% hObject    handle to hold_axis_lims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold_axis_lims

% Define the event function
function [VALUE, ISTERMINAL, DIRECTION] = MyEventFunction(T, Y)
%The event function stops the intergration is VALUE == 0 and 
%ISTERMINAL==1

%a. Define the timeout in seconds
TimeOut = 0.5;
%
%b. The solver runs until this VALUE is negative (does not change the sign)
    VALUE = toc-TimeOut;


%c. The function should terminate the execution, so
ISTERMINAL = 1;

%d. The direction does not matter
DIRECTION = 0;

% what is funny, it works!! 


% --- Executes on button press in IC_default.
function IC_default_Callback(hObject, eventdata, handles)
% hObject    handle to IC_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.initialConditions = [0.5 0.5];
% set(handles.initCondDisp, 'String', 'Initial Conditions: [0.5 0.5]');
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over IC_default.
function IC_default_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to IC_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in timespan_default.
function timespan_default_Callback(hObject, eventdata, handles)
% hObject    handle to timespan_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.timeSpan = 50;
% set(handles.timeSpanDisp, 'String', 'Time Span: 50s');
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over timespan_default.
function timespan_default_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to timespan_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function Theta_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% coords = get(hObject, 'CurrentPoint');
% coords = min(max(coords,0),1);
% handles.initialConditions = [coords(1,1), coords(1,2)];
% ic = handles.initialConditions;
% icstr = sprintf('Initial Conditions: [%.3g %.3g]', ic(1), ic(2));
% set(handles.initCondDisp, 'String', icstr);
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)


% Attempts to measure the period of a limit cycle
function period = measureLimitCycle(time, biomass)
% peakLocs = peakfinder(biomass);
% peakTimes = time(peakLocs);
% period = mean(diff(peakTimes));






% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Run.
function Run_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thetar = 1; % radians
thetardot = 0;
zeta = 1;
omega = 10;
Jhat = 1;
J = 2;
mup = 10;
mud = 10;

t = 0:0.01:100;
[tout, ~, yout] = sim('satelliteProject', t);
goal = tout./tout.*thetar;
x1 = yout(:,5); 
x2 = yout(:,6); 
u = yout(:,7); 
input = yout(:,8); 

% figure(2)
% clf
% plot(tout,input);

% %Plot time domain results
axes(handles.Dhat)

% % Store old x and y axis limits
% old_xlims = xlim;
% old_ylims = ylim;

cla(handles.Input);
% hold all
plot(tout,input);
% plot(t,x(:,1), 'green');
% plot(t,x(:,2), 'red');
% legend('Prey','Predator');

T = 0.01; %s
I = T*sum(x1.^2+x2.^2+(u.^2)/100);

% disp('I: ')
% disp(I);

% figure(1)
% clf
% subplot(3,1,1);
% plot(tout,yout(:,1)); % Theta
% hold all
% plot(tout,yout(:,4)); % Theta_m
% % plot(tout,goal,'--'); % Theta_r
% % ylim([0.04 0.07]);
% xlabel('Time (s)');
% ylabel('Theta (rad)');
% legend('Theta','Theta_m','Theta_r');
% title('Theta over time');
% subplot(3,1,2);
% plot(tout,yout(:,2)); % P
% hold all
% % plot(tout,tout./tout*J*omega^2); % P_hat
% xlabel('Time (s)');
% ylabel('P');
% title('P over time');
% subplot(3,1,3); % D
% plot(tout,yout(:,3));
% hold all
% % plot(tout,tout./tout*J*2*zeta*omega); % D_hat
% xlabel('Time (s)');
% ylabel('D');
% title('D over time');











