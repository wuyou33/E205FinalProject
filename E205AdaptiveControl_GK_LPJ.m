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

% Last Modified by GUIDE v2.5 16-Dec-2014 23:16:55

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
handles.initialAngle = 0.5;
handles.initialVel = 0.5;
handles.muP = 10;
handles.muD = 10;
handles.ditherAmplitude = 0.01;
handles.ditherFreq = 100;
handles.timeSpan = 50;
handles.plotSpan = 10;
handles.refSignal = 1;
handles.refAmplitude = 1;
handles.refFreq = 100;
handles.results = {};

set(handles.sliderDitherAmp,'value', log10(handles.ditherAmplitude));
set(handles.sliderDitherFreq,'value', log10(handles.ditherFreq));

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
% set(handles.sliderDitherAmp,'value', mu);
% set(handles.sliderDitherFreq,'value', sigma);
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
% hObject    handle to stabilityBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Updates values of a and b (reaction constants)

function sliderDitherAmp_Callback(hObject, eventdata, handles)
% hObject    handle to dispDitherAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispDitherAmp as text
%        str2double(get(hObject,'String')) returns contents of dispDitherAmp as a double
sliderValue = get(hObject,'value');
handles.ditherAmplitude = 10^sliderValue;
ditherAmp_str = sprintf('%g', handles.ditherAmplitude);
set(handles.dispDitherAmp, 'String', ditherAmp_str);
% handles = updatePhasePlot(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function sliderDitherAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispDitherAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sliderDitherFreq_Callback(hObject, eventdata, handles)
% hObject    handle to dispDitherFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispDitherFreq as text
%        str2double(get(hObject,'String')) returns contents of dispDitherFreq as a double
% handles.sigma = get(hObject,'value');
% sigma_str = sprintf('%g', handles.sigma);
% set(handles.dispDitherFreq, 'String', sigma_str);
% handles = updatePhasePlot(handles);
% guidata(hObject, handles)
sliderValue = get(hObject,'value');
handles.ditherFreq = 10^sliderValue;
ditherFreq_str = sprintf('%g', handles.ditherFreq);
set(handles.dispDitherFreq, 'String', ditherFreq_str);
% handles = updatePhasePlot(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function sliderDitherFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispDitherFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dispDitherAmp_Callback(hObject, eventdata, handles)
% hObject    handle to dispDitherAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispDitherAmp as text
%        str2double(get(hObject,'String')) returns contents of dispDitherAmp as a double
newditherAmp = str2double(get(hObject,'String'));
minditherAmp = 0;
maxditherAmp = 1;
if isnan(newditherAmp)
    set(hObject, 'String', '0.01')
    newditherAmp = 0.01;
elseif newditherAmp > maxditherAmp
    set(hObject, 'String', num2str(maxditherAmp))
    newditherAmp = maxditherAmp;
elseif newditherAmp < minditherAmp
    set(hObject, 'String', num2str(minditherAmp))
    newditherAmp = minditherAmp;
end
handles.ditherAmplitude = newditherAmp;
set(handles.sliderDitherAmp,'value', log10(handles.ditherAmplitude));
% handles = updatePhasePlot(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispDitherAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispDitherAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dispDitherFreq_Callback(hObject, eventdata, handles)
% hObject    handle to dispDitherFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispDitherFreq as text
%        str2double(get(hObject,'String')) returns contents of dispDitherFreq as a double
newditherFreq = str2double(get(hObject,'String'));

minditherFreq = 1;
maxditherFreq = 1000;
if isnan(newditherFreq)
    set(hObject, 'String', '0.01')
    newditherFreq = 0.01;
elseif newditherFreq > maxditherFreq
    set(hObject, 'String', num2str(maxditherFreq))
    newditherFreq = maxditherFreq;
elseif newditherFreq < minditherFreq
    set(hObject, 'String', num2str(minditherFreq))
    newditherFreq = minditherFreq;
end
handles.ditherFreq = newditherFreq;
set(handles.sliderDitherFreq,'value', log10(handles.ditherFreq));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function dispDitherFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispDitherFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function dispInitAngle_Callback(hObject, eventdata, handles)
% hObject    handle to dispInitAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispInitAngle as text
%        str2double(get(hObject,'String')) returns contents of dispInitAngle as a double
newAngle = str2double(get(hObject,'String'));

minAngle = 0;
maxAngle = 2*pi;
if isnan(newAngle)
    set(hObject, 'String', '0.5')
    newAngle = 0.5;
elseif newAngle > maxAngle
    set(hObject, 'String', num2str(maxAngle))
    newAngle = maxAngle;
elseif newAngle < minAngle
    set(hObject, 'String', num2str(minAngle))
    newAngle = minAngle;
end
handles.initialAngle = newAngle;
% handles = updatePhasePlot(handles);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function dispInitAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispInitAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dispTimeSpan_Callback(hObject, eventdata, handles)
% hObject    handle to dispTimeSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispTimeSpan as text
%        str2double(get(hObject,'String')) returns contents of dispTimeSpan as a double
newSpan = str2double(get(hObject,'String'));

minSpan = 1;
maxSpan = 200;
if isnan(newSpan)
    set(hObject, 'String', '50')
    newSpan = 50;
elseif newSpan > maxSpan
    set(hObject, 'String', num2str(maxSpan))
    newSpan = maxSpan;
elseif newSpan < minSpan
    set(hObject, 'String', num2str(minSpan))
    newSpan = minSpan;
end
handles.timeSpan = newSpan;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispTimeSpan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispTimeSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on button press in defaultInitAngle.
function defaultInitAngle_Callback(hObject, eventdata, handles)
% hObject    handle to defaultInitAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.initialAngle = 0.5;
set(handles.dispInitAngle, 'String', '0.5');
guidata(hObject, handles)

% --- Executes on button press in defaultTimeSpan.
function defaultTimeSpan_Callback(hObject, eventdata, handles)
% hObject    handle to defaultTimeSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.timeSpan = 50;
set(handles.dispTimeSpan, 'String', '50');
guidata(hObject, handles)

function dispInertia_Callback(hObject, eventdata, handles)
% hObject    handle to dispInertia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispInertia as text
%        str2double(get(hObject,'String')) returns contents of dispInertia as a double
newInertia = str2double(get(hObject,'String'));

minInertia = 0.1;
maxInertia = 100;
if isnan(newInertia)
    set(hObject, 'String', '2')
    newInertia = 2;
elseif newInertia > maxInertia
    set(hObject, 'String', num2str(maxInertia))
    newInertia = maxInertia;
elseif newInertia < minInertia
    set(hObject, 'String', num2str(minInertia))
    newInertia = minInertia;
end
handles.Jactual = newInertia;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispInertia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispInertia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultInertia.
function defaultInertia_Callback(hObject, eventdata, handles)
% hObject    handle to defaultInertia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Jactual = 2;
set(handles.dispInertia, 'String', '2');
guidata(hObject, handles)


function dispInertiaGuess_Callback(hObject, eventdata, handles)
% hObject    handle to dispInertiaGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispInertiaGuess as text
%        str2double(get(hObject,'String')) returns contents of dispInertiaGuess as a double
newInertia = str2double(get(hObject,'String'));

minInertia = 0.1;
maxInertia = 100;
if isnan(newInertia)
    set(hObject, 'String', '1')
    newInertia = 2;
elseif newInertia > maxInertia
    set(hObject, 'String', num2str(maxInertia))
    newInertia = maxInertia;
elseif newInertia < minInertia
    set(hObject, 'String', num2str(minInertia))
    newInertia = minInertia;
end
handles.Jguess = newInertia;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispInertiaGuess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispInertiaGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultInertiaGuess.
function defaultInertiaGuess_Callback(hObject, eventdata, handles)
% hObject    handle to defaultInertiaGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Jguess = 1;
set(handles.dispInertiaGuess, 'String', '1');
guidata(hObject, handles)


function dispPRate_Callback(hObject, eventdata, handles)
% hObject    handle to dispPRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispPRate as text
%        str2double(get(hObject,'String')) returns contents of dispPRate as a double
newRate = str2double(get(hObject,'String'));

minRate = 0.1;
maxRate = 100;
if isnan(newRate)
    set(hObject, 'String', '10')
    newRate = 2;
elseif newRate > maxRate
    set(hObject, 'String', num2str(maxRate))
    newRate = maxRate;
elseif newRate < minRate
    set(hObject, 'String', num2str(minRate))
    newRate = minRate;
end
handles.muP = newRate;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispPRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispPRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultPRate.
function defaultPRate_Callback(hObject, eventdata, handles)
% hObject    handle to defaultPRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.muP = 10;
set(handles.dispPRate, 'String', '10');
guidata(hObject, handles)


function dispDRate_Callback(hObject, eventdata, handles)
% hObject    handle to dispDRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispDRate as text
%        str2double(get(hObject,'String')) returns contents of dispDRate as a double
newRate = str2double(get(hObject,'String'));

minRate = 0.1;
maxRate = 100;
if isnan(newRate)
    set(hObject, 'String', '10')
    newRate = 2;
elseif newRate > maxRate
    set(hObject, 'String', num2str(maxRate))
    newRate = maxRate;
elseif newRate < minRate
    set(hObject, 'String', num2str(minRate))
    newRate = minRate;
end
handles.muD = newRate;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispDRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispDRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultDRate.
function defaultDRate_Callback(hObject, eventdata, handles)
% hObject    handle to defaultDRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.muD = 10;
set(handles.dispDRate, 'String', '10');
guidata(hObject, handles)


function dispInitVel_Callback(hObject, eventdata, handles)
% hObject    handle to dispInitVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispInitVel as text
%        str2double(get(hObject,'String')) returns contents of dispInitVel as a double
newVel = str2double(get(hObject,'String'));

minVel = -1000;
maxVel = 1000;
if isnan(newVel)
    set(hObject, 'String', '0.5')
    newVel = 2;
elseif newVel > maxVel
    set(hObject, 'String', num2str(maxVel))
    newVel = maxVel;
elseif newVel < minVel
    set(hObject, 'String', num2str(minVel))
    newVel = minVel;
end
handles.initialVel = newVel;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispInitVel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispInitVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultInitVel.
function defaultInitVel_Callback(hObject, eventdata, handles)
% hObject    handle to defaultInitVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.initialVel = 0.5;
set(handles.dispInitVel, 'String', '0.5');
guidata(hObject, handles)


% --- Executes on mouse press over axes background.
function Theta_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make these inputs?
zeta = 1;
omega = 10;

assignin('base', 'Amplitude', handles.refAmplitude);
assignin('base', 'Frequency', handles.refFreq);
assignin('base', 'Period', 5);
assignin('base', 'zeta', zeta);
assignin('base', 'omega', omega);
assignin('base', 'Jhat', handles.Jguess);
assignin('base', 'J', handles.Jactual);
assignin('base', 'mup', handles.muP);
assignin('base', 'mud', handles.muD);
assignin('base', 'inputChoice', handles.refSignal);
assignin('base', 'ditherAmp', handles.ditherAmplitude);
assignin('base', 'ditherFreq', handles.ditherFreq);

t = 0:0.01:handles.timeSpan;
[tout, ~, yout] = sim('satelliteProject', t);

x1 = yout(:,5); 
x2 = yout(:,6); 
u = yout(:,7); 
input = yout(:,8); 

axes(handles.InputSignal);
cla(handles.InputSignal);
plot(tout,input);

axes(handles.Theta);
cla(handles.Theta);
plot(tout,yout(:,1)); % Theta
hold all
plot(tout,yout(:,4)); % Theta_m
xlabel('Time (s)');
ylabel('\theta (rad)');
legend('\theta','\theta_m');
title('Theta over time');

axes(handles.Phat);
cla(handles.Phat);
plot(tout,yout(:,2)); % P
hold all
P_hat = tout./tout*handles.Jactual*omega^2;
plot(tout,P_hat, 'r');
ylim([min(yout(:,2)), max([yout(:,2);P_hat])*1.5]);
xlabel('Time (s)');
ylabel('P');
legend('P','P\_hat');
title('P over time');

axes(handles.Dhat);
cla(handles.Dhat);
plot(tout,yout(:,3));
hold all
D_hat = tout./tout*handles.Jactual*2*zeta*omega;
plot(tout, D_hat, 'r'); % D_hat
ylim([min(yout(:,3)), max([yout(:,3);D_hat])*1.5]);
xlabel('Time (s)');
ylabel('D');
legend('D','D\_hat');
title('D over time');


guidata(hObject, handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Run.
function Run_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function dispSignalAmp_Callback(hObject, eventdata, handles)
% hObject    handle to dispSignalAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispSignalAmp as text
%        str2double(get(hObject,'String')) returns contents of dispSignalAmp as a double
newSignalAmp = str2double(get(hObject,'String'));

minSignalAmp = 1;
maxSignalAmp = 2*pi;
if isnan(newSignalAmp)
    set(hObject, 'String', '1')
    newSignalAmp = 1;
elseif newSignalAmp > maxSignalAmp
    set(hObject, 'String', num2str(maxSignalAmp))
    newSignalAmp = maxSignalAmp;
elseif newSignalAmp < minSignalAmp
    set(hObject, 'String', num2str(minSignalAmp))
    newSignalAmp = minSignalAmp;
end
handles.refAmplitude = newSignalAmp;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispSignalAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispSignalAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dispSignalFreq_Callback(hObject, eventdata, handles)
% hObject    handle to dispSignalFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dispSignalFreq as text
%        str2double(get(hObject,'String')) returns contents of dispSignalFreq as a double
newSignalFreq = str2double(get(hObject,'String'));

minSignalFreq = 1;
maxSignalFreq = 1000;
if isnan(newSignalFreq)
    set(hObject, 'String', '100')
    newSignalFreq = 100;
elseif newSignalFreq > maxSignalFreq
    set(hObject, 'String', num2str(maxSignalFreq))
    newSignalFreq = maxSignalFreq;
elseif newSignalFreq < minSignalFreq
    set(hObject, 'String', num2str(minSignalFreq))
    newSignalFreq = minSignalFreq;
end
handles.refFreq = newSignalFreq;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function dispSignalFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispSignalFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in signalTypeSelect.
function signalTypeSelect_Callback(hObject, eventdata, handles)
% hObject    handle to signalTypeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns signalTypeSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from signalTypeSelect
contents = cellstr(get(hObject,'String'));
selected = contents{get(hObject,'Value')};
if (strcmp(selected, 'Constant (Set Point)'))
    set(handles.textSignalAmp, 'String', 'Set Point (rad)')
    set(handles.textSignalFreq, 'Visible', 'off')
    set(handles.dispSignalFreq, 'Visible', 'off')
    set(handles.defaultSignalFreq, 'Visible', 'off')
    handles.refSignal = 1;
elseif (strcmp(selected, 'Constant Rotation (Ramp)'))
    set(handles.textSignalAmp, 'String', 'Set Point (rad)')
    set(handles.textSignalFreq, 'Visible', 'off')
    set(handles.dispSignalFreq, 'Visible', 'off')
    set(handles.defaultSignalFreq, 'Visible', 'off')
    handles.refSignal = 4;
elseif (strcmp(selected, 'Sinusoid'))
    set(handles.textSignalAmp, 'String', 'Amplitude (rad)')
    set(handles.textSignalFreq, 'Visible', 'on')
    set(handles.dispSignalFreq, 'Visible', 'on')
    set(handles.defaultSignalFreq, 'Visible', 'on')
    handles.refSignal = 3;
else %(strcmp(selected, 'Pulse Train'))
    set(handles.textSignalAmp, 'String', 'Amplitude (rad)')
    set(handles.textSignalFreq, 'Visible', 'on')
    set(handles.dispSignalFreq, 'Visible', 'on')
    set(handles.defaultSignalFreq, 'Visible', 'on')
    handles.refSignal = 2;
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function signalTypeSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to signalTypeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultSignalAmplitude.
function defaultSignalAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to defaultSignalAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.refAmplitude = 1;
ampstr = sprintf('%g', handles.refAmplitude);
set(handles.dispSignalAmp, 'String', ampstr);
guidata(hObject, handles)

% --- Executes on button press in defaultSignalFreq.
function defaultSignalFreq_Callback(hObject, eventdata, handles)
% hObject    handle to defaultSignalFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.refFreq = 100;
freqstr = sprintf('%g', handles.refFreq);
set(handles.dispSignalFreq, 'String', freqstr);
guidata(hObject, handles)


