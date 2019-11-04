
function varargout = nonrigid_assited_image_distortion(varargin)
	% NONRIGID_ASSITED_IMAGE_DISTORTION MATLAB code for nonrigid_assited_image_distortion.fig
	%      NONRIGID_ASSITED_IMAGE_DISTORTION, by itself, creates a new NONRIGID_ASSITED_IMAGE_DISTORTION or raises the existing
	%      singleton*.
	%
	%      H = NONRIGID_ASSITED_IMAGE_DISTORTION returns the handle to a new NONRIGID_ASSITED_IMAGE_DISTORTION or the handle to
	%      the existing singleton*.
	%
	%      NONRIGID_ASSITED_IMAGE_DISTORTION('CALLBACK',hObject,eventData,handles,...) calls the local
	%      function named CALLBACK in NONRIGID_ASSITED_IMAGE_DISTORTION.M with the given input arguments.
	%
	%      NONRIGID_ASSITED_IMAGE_DISTORTION('Property','Value',...) creates a new NONRIGID_ASSITED_IMAGE_DISTORTION or raises the
	%      existing singleton*.  Starting from the left, property value pairs are
	%      applied to the GUI before nonrigid_assited_image_distortion_OpeningFcn gets called.  An
	%      unrecognized property name or invalid value makes property application
	%      stop.  All inputs are passed to nonrigid_assited_image_distortion_OpeningFcn via varargin.
	%
	%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
	%      instance to run (singleton)".
	%
	% See also: GUIDE, GUIDATA, GUIHANDLES

	% Edit the above text to modify the response to help nonrigid_assited_image_distortion

	% Last Modified by GUIDE v2.5 06-Jun-2019 21:22:36

	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 0;
	gui_State = struct('gui_Name',       mfilename, ...
	                   'gui_Singleton',  gui_Singleton, ...
	                   'gui_OpeningFcn', @nonrigid_assited_image_distortion_OpeningFcn, ...
	                   'gui_OutputFcn',  @nonrigid_assited_image_distortion_OutputFcn, ...
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

% --- Executes just before nonrigid_assited_image_distortion is made visible.
function nonrigid_assited_image_distortion_OpeningFcn(hObject, eventdata, handles, varargin)
	% This function has no output args, see OutputFcn.
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	% varargin   command line arguments to nonrigid_assited_image_distortion (see VARARGIN)

	% Choose default command line output for nonrigid_assited_image_distortion
	handles.output = hObject;

	% Update handles structure
	guidata(hObject, handles);

	% UIWAIT makes nonrigid_assited_image_distortion wait for user response (see UIRESUME)
	% uiwait(handles.figure1);
	try
	    functionname='nonrigid_assited_image_distortion.m';
	    functiondir=which(functionname);
	    functiondir=functiondir(1:end-length(functionname));
	    addpath([functiondir '/functions'])
	    addpath([functiondir '/functions_nonrigid'])
	catch me
	    disp(me.message);
	end

	% Disable warning
	warning('off', 'MATLAB:maxNumCompThreads:Deprecated')
	set(handles.axes1, 'xtick', [], 'ytick', [], 'Box', 'on');
	%set(handles.axes2, 'xtick', [], 'ytick', [], 'Box', 'on');
	%axis(handles.axes1,'off');
	%axis(handles.axes2,'off');

	%--- Initializations
	params.empty = true;
	params.option = 'build_grid';
	params.handles = handles;
	params.optvideo.status = false;
	params.optvideo.pause = false;
	params.alpha = 0.5;
	params.padding = 2;
	params.channel = 1; % red channel by default
	set_params(params);
	set_menu_checks();

% --- Outputs from this function are returned to the command line.
function varargout = nonrigid_assited_image_distortion_OutputFcn(hObject, eventdata, handles) 
	% varargout  cell array for returning output args (see VARARGOUT);
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;

% --------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
	% hObject    handle to menu_file (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------
function open_image_Callback(hObject, eventdata, handles)
	% hObject    handle to open_image (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	disp('open_image');
	params = get_params();

	[file, fpath] = uigetfile( ...
		{'*.jpg';'*.jepg';'*.png';'*.bmp';'*.tiff';'*.*'}, ...
		'Open the image' );
	if isequal(file,0)
		disp('User selected Cancel')
		return
	else
		im = imread(fullfile(fpath, file));
	end
	%set_params(params);
	%try
	if strcmpi(params.option, 'build_grid')
		initialization(im);
	elseif strcmpi(params.option, 'apply_grid')
		initialization2(im);
	end
	%catch
	%	return
	%end

function set_params(params)
	%--- Set params struct in figure
	setappdata(gcf, 'paramstransformix', params);

function params = get_params()
	%--- Get params struct stored in figure
	params = getappdata(gcf, 'paramstransformix');

function initialization(im)
	params = get_params();
	if isfield(params, 'gui_timer')
		stop(params.gui_timer); 
		delete(params.gui_timer);
		try
			close(params.f);
		catch
			disp('Second figure already was closed');
		end
	end

	params.empty = false;
	params.optvideo.status = false;
	params.optvideo.pause = true;
	params.im = im2double(im);
	imsize = size(params.im);
	imsizepadding = imsize;

	% Get the power of max refinement steps
	powermax = min( floor( log2(imsize(1:2)*.25) ) );

	% Get b-spline grid spacing in x and y direction
	params.spacing = [2^powermax, 2^powermax];
	params.uniform = make_init_grid( params.spacing, ...
		imsize(1:2) );
	params.grid = params.uniform;

	if length(imsize) == 3
		params.isrgb = true;
		params.imrgb = params.im;

		params.imred = params.im;
		params.imred(:,:,2:3) = 0.0;
		params.imblue = params.im;
		params.imblue(:,:,1:2) = 0.0;
		params.imgreen = params.im;
		params.imgreen(:,:,1:2:3) = 0.0;

		params.uniformrgb = params.uniform;
		params.uniformred = params.uniform;
		params.uniformblue = params.uniform;
		params.uniformgreen = params.uniform;

		params.gridrgb = params.grid;
		params.gridred = params.grid;
		params.gridblue = params.grid;
		params.gridgreen = params.grid;

		params.spacingrgb = params.spacing;
		params.spacingred = params.spacing;
		params.spacingblue = params.spacing;
		params.spacinggreen = params.spacing;

		switch(params.channel)
			case 0
		        params.im = params.imrgb;
		    case 1
		        params.im = params.imred;
		    case 2
		        params.im = params.imgreen;
		    case 3
		        params.im = params.imblue;
		    otherwise 
		end
	else
		params.isrgb = false;
	end

	params.handle_grid_axes = params.handles.axes1;
	params.f = figure('units','normalized','outerposition',[0 0 1 1]);
	params.h1 = axes(params.f, 'Position', [.05 .1 .425 .8]);
	params.h2 = axes(params.f, 'Position', [.525 .1 .425 .8]);
	params.handle_control_point_selected = [];
	params.view_option = 1;
	params.view_swap = 0;

	imsizepadding(1:2) = imsize(1:2) + params.padding * 2;
	params.imt = params.im;	
	params.ims = zeros(imsizepadding);
	params.imsgrid = zeros(imsizepadding);

	x1 = params.padding + 1; 
	x2 = params.padding + imsize(1);
	y1 = params.padding + 1; 
	y2 = params.padding + imsize(2);

	params.ims( x1:x2, y1:y2, : ) = params.imt;
	params.imsgrid( x1:x2, y1:y2, : ) = params.imt;
	
	igrid = params.grid + params.padding;
	params.imsgrid = get_checkerboard(igrid, params.imsgrid, params.alpha);

	%params.handle_imgrid = imshow(params.imsgrid, 'Parent', ...
	%	params.handle_grid_axes);
	%size(params.imsgrid)
	%size(params.imch)
	colormap(params.handle_grid_axes, gray)
	params.handle_imgrid = image(params.imsgrid, 'Parent', ...
		params.handle_grid_axes, 'CDataMapping','scaled');
	set(params.handle_grid_axes, 'xtick', [], 'ytick', []);

	params.handle_ims = imshow(params.ims, 'Parent', ...
		params.h1);
	params.handle_imt = imshow(params.imt, 'Parent', ...
		params.h2);

	params.handle_control_points = [];
	%params.handles = handles;

	figure(params.handles.figure1)
	params.gui_timer = timer('TimerFcn', ...
		'nonrigid_assited_image_distortion(''TimerFcn'',gcf,[],guidata(gcf));','Period',1,'ExecutionMode','fixedDelay');
	set_params(params);

	start(params.gui_timer);
	% ToDo: setting Jstatic
	bspline_tranformix_image();
	show_slices();
	show_grid();

	figure(params.f)
	set_params(params);
	figure(params.handles.figure1)

function initialization2(im)
	params = get_params();
	if isfield(params, 'gui_timer')
		stop(params.gui_timer); 
		delete(params.gui_timer);
		try
			close(params.f);
		catch
			disp('Second figure already was closed');
		end
	end

	% Get the power of max refinement steps
	imsize = size(im);
	powermax = min( floor( log2(imsize(1:2)*.25) ) );

	% Get b-spline grid spacing in x and y direction
	params.spacing = [2^powermax, 2^powermax];
	params.uniform = make_init_grid( params.spacing, ...
		imsize(1:2) );
	params.grid = params.uniform;

	params.empty = false;
	params.optvideo.status = false;
	params.optvideo.pause = true;

	if length(imsize) == 3
		params.isrgb = true;
		params.imrgb = im2double(im);
		params.im = params.imrgb;
		switch(params.channel)
			case 0
		        params.im = params.imrgb;
		    case 1
				params.im(:,:,2:3) = 0.0;
		    case 2
		        params.im(:,:,1:2:3) = 0.0;
		    case 3
		        params.im(:,:,1:2) = 0.0;
		    otherwise 
		end
	else
		params.isrgb = false;
		params.im = im2double(im);
	end
		
	
	params.handle_grid_axes = params.handles.axes1;
	params.imt = [];
	params.handle_imgrid = imshow(params.im, 'Parent', ...
		params.handle_grid_axes);
	%params.handle_imgrid = image(params.im, 'Parent', ...
	%	params.handle_grid_axes);
	params.handle_control_points = [];
	params.handle_control_point_selected = [];

	set_params(params);
	%bspline_tranformix_image();

function bspline_tranformix_image()
	disp('bspline_tranformix_image');
	params = get_params();
	params.real_grid = params.uniform + ( ...
		params.uniform - params.grid );
	params.imt = bspline_transform( ...
		params.real_grid, ...
		params.im, ...
		params.spacing, ...
		3 );

	if ~strcmpi(params.option, 'apply_grid')	
		switch(params.channel)
			case 0
		        params.imrgb = params.im;
		        params.uniformrgb = params.uniform;
		        params.gridrgb = params.grid;
		        params.spacingrgb = params.spacing;
		    case 1
		        params.imred = params.im;
		        params.uniformred = params.uniform;
		        params.gridred = params.grid;
		        params.spacingred = params.spacing;
		    case 2
		        params.imgreen = params.im;
		        params.uniformgreen = params.uniform;
		        params.gridgreen = params.grid;
		        params.spacinggreen = params.spacing;
		    case 3
		        params.imblue = params.im;
		        params.uniformblue = params.uniform;
		        params.gridblue = params.grid;
		        params.spacingblue = params.spacing;
		    otherwise 
		end
	end

	set_params(params);

	if ~strcmpi(params.option, 'apply_grid')
		figure(params.f)
		set_params(params);
		figure(params.handles.figure1)
	end

function bspline_tranformix_frame()
	%disp('bspline_tranformix_frames');
	params = get_params();
	params.real_grid = params.uniform + ( ...
		params.uniform - params.grid );
	params.framet = bspline_transform( ...
		params.real_grid, ...
		params.frame, ...
		params.spacing, ...
		3 );
	set_params(params);

function show_video()
	params = get_params();

	while hasFrame(params.video)
		params = get_params();
		if isempty(params), return, end
		if params.optvideo.pause
	    	break
	    end

	    params.frame = im2double(readFrame(params.video));
	    switch(params.channel)
			case 0
		        params.frame = params.frame;
		    case 1
				params.frame(:,:,2:3) = 0.0;
		    case 2
		        params.frame(:,:,1:2:3) = 0.0;
		    case 3
		        params.frame(:,:,1:2) = 0.0;
		    otherwise 
		end

		set_params(params);

		bspline_tranformix_frame()
	    set(params.handle_imgrid, 'CData', params.framet);
	    pause(.5*(1/params.video.FrameRate));
	    
	end

function show_grid()
	params = get_params();

	if(ishandle(params.handle_control_points))
	    delete(params.handle_control_points);
	    % delete(params.handle_grid_vlines);
	    % delete(params.handle_grid_hlines);
	end
	gridsize = size(params.grid);
	numelgrid = prod(gridsize(1:2));
	params.handle_control_points = zeros(numelgrid, 1);
	params.handle_grid_vlines = zeros(numelgrid, 1);
	params.handle_grid_hlines = zeros(numelgrid, 1);

	hold(params.handle_grid_axes,'on');
	igrid = params.grid + params.padding;
	k = 0;
	for i=1:gridsize(1),
		for j=1:gridsize(2),
		    k=k+1;
		    hmin = min(j+1,size(igrid,2));
		    vmin = min(i+1,size(igrid,1));
		    % params.handle_grid_hlines(k) = plot( ...
		    % 	params.handle_grid_axes, ...
		    % 	[igrid(i,j,2) igrid(i,hmin,2)],...
		    % 	[igrid(i,j,1) igrid(i,hmin,1)], '--');
		    % params.handle_grid_vlines(k) = plot( ...
		    % 	params.handle_grid_axes, ...
		    % 	[igrid(i,j,2) igrid(vmin,j,2)], ...
		    % 	[igrid(i,j,1) igrid(vmin,j,1)], '--');
		    params.handle_control_points(k) = plot( ...
		    	params.handle_grid_axes, ...
		    	igrid(i,j,2), igrid(i,j,1), ...
		    	'wo','MarkerSize',8);%, 'MarkerFaceColor', 'w');

		    % if (j == 1) || (j == gridsize(2) - 2)
		    % 	set(params.handle_grid_hlines(k), ...
			   %  	'XData', [], ...
			   %  	'YData', [] )
		    % end
		    % if (i == 1) || (i == gridsize(1) - 2)
		    % 	set(params.handle_grid_vlines(k), ...
			   %  	'XData', [], ...
			   %  	'YData', [] )
		    % end
		end
	end

	set(params.handle_control_points(:),'ButtonDownFcn','nonrigid_assited_image_distortion(''pointGridButtonDownFcn'',gcbo,[],guidata(gcbo))');

	hold(params.handle_grid_axes,'off');
	set_params(params);

function show_slices()
	disp('show_slices');
	params = get_params();
	try
		imsize = size(params.imt);
	catch me
	    return
	end

	x1 = params.padding + 1; 
	x2 = params.padding + imsize(1);
	y1 = params.padding + 1; 
	y2 = params.padding + imsize(2);

	if(params.view_option == 1)
		if(params.view_swap == 0)
			params.ims( x1:x2, y1:y2, : ) = .85*params.imt;
		else
			params.ims( x1:x2, y1:y2, : ) = params.im;
		end
	end
	params.imsgrid( x1:x2, y1:y2, : ) = params.imt;

	igrid = params.grid + params.padding;
	params.imsgrid = get_checkerboard(igrid, params.imsgrid, params.alpha);

	%params.ims(params.ims<0)=0;
	%params.ims(params.ims>1)=1;

	set(params.handle_imgrid, 'CData', params.imsgrid);
	set(params.handle_ims, 'CData', params.ims);
	set(params.handle_imt, 'CData', params.imt);
	%drawnow

function [imch] = get_checkerboard(igrid, im, alphaa)
	gridsize = size(igrid);
	imsize = size(im);
	imch = im;
	%imch = zeros(imsize(1:2));
	colour = {'black', 'white'};
	for i1 = 2:gridsize(1) - 3
		val = mod(i1, 2);
		for i2 = 2:gridsize(2) - 3
			x1 = round(igrid(i1,i2,2));
			y1 = round(igrid(i1,i2,1));
			x2 = round(min(igrid(i1,i2+1,2), imsize(2)));
			y2 = round(min(igrid(i1,i2+1,1), imsize(2)));
			x3 = round(min(igrid(i1+1,i2+1,2), imsize(2)));
			y3 = round(min(igrid(i1+1,i2+1,1), imsize(1)));
			x4 = round(min(igrid(i1+1,i2,2), imsize(2)));
			y4 = round(min(igrid(i1+1,i2,1), imsize(1)));
			imch = insertShape(imch,'FilledPolygon',...
        		[x1,y1, x2,y2, x3,y3, x4,y4],...
        		'Color', colour{val+1}, 'Opacity', alphaa);
			%imch(r1:r2, c1:c2) = val;
			if val == 1
				val = 0;
			else
				val = 1;
			end
		end
	end

	if length(imsize) < 3
		imch = imch(:,:,1);
	end

	% if length(imsize) == 3
	% 	imch3 = zeros(imsize);
	% 	for i3 = 1:3
	% 		imch3(:,:,i3) = imch;
	% 	end
	% 	imch = imch3;
	% end

function cursor_position_update()
	%disp('cursor_position_update');
	params = get_params();
	if( isempty(params) ), return, end
	if( params.empty ), return, end
	p = get(params.handle_grid_axes, 'CurrentPoint');
	params.mouse_position=[p(1, 1) p(1, 2)];
	set_params(params);

function update_grid()
	params = get_params();
	gridsize = size(params.grid);
	[j,i] = ind2sub( [gridsize(2), gridsize(1)], ...
		params.control_point_selected );
	params.grid(i,j,1) = params.mouse_position(2) - params.padding;
	params.grid(i,j,2) = params.mouse_position(1) - params.padding;
	params.handle_control_point_selected = [];
	set_params(params);
	bspline_tranformix_image();
	show_slices();

	%show_grid();
	%k = params.control_point_selected;
    %set(params.handle_control_points(k), ...
    %	'XData', params.mouse_position(1));
    %set(params.handle_control_points(k), ...
    %	'YData', params.mouse_position(2));	

function set_menu_checks()
	params = get_params();
	if strcmpi(params.option, 'build_grid')
		set(params.handles.build_grid, 'Checked', 'on');
		set(params.handles.apply_grid, 'Checked', 'off');
	elseif strcmpi(params.option, 'apply_grid')
		set(params.handles.build_grid, 'Checked', 'off');
		set(params.handles.apply_grid, 'Checked', 'on');
	else
		disp('Another condition');
	end

function clear_figure1()
	params = get_params();

	cla(params.handles.axes1);
	axis(params.handles.axes1,'on');
	set(params.handles.axes1, 'xtick', [], 'ytick', [], 'Box', 'on');

	if isfield(params, 'gui_timer')
		stop(params.gui_timer);
		delete(params.gui_timer);
		try
			close(params.f);
		catch
			disp('Second figure already was closed');
		end
	end
		

function TimerFcn(hObject, eventdata, handles)
	disp('TimerFcn');
	params = get_params();
	if isempty(params), return, end
	if(params.empty), return, end
	if(params.view_option == 1) 
		params.view_swap = params.view_swap + 1;
		if( params.view_swap > 1 ), params.view_swap = 0; end
		set_params(params);
		show_slices();
	end

function pointGridButtonDownFcn(hObject, eventdata, handles)
	disp('pointGridButtonDownFcn');
	params = get_params(); 
	if( isempty(params) ), return, end
	if(params.empty), return, end
	if (strcmpi(params.option, 'apply_grid')), return, end
	params.control_point_selected = find( ...
		params.handle_control_points == gcbo);
	params.handle_control_point_selected = gcbo;
	set(params.handle_control_point_selected, ...
		'MarkerSize', 10, 'MarkerFaceColor', 'k');
	set_params(params);

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
	% hObject    handle to figure1 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	%disp('figure1_WindowButtonMotionFcn');
	cursor_position_update();
	params = get_params();
	if( isempty(params) ), return, end
	if (params.empty), return, end
	if (strcmpi(params.option, 'apply_grid')), return, end
	if(~isempty(params.handle_control_point_selected))
		k = params.control_point_selected;
	    set(params.handle_control_point_selected, ...
	    	'XData', params.mouse_position(1));
	    set(params.handle_control_point_selected, ...
	    	'YData', params.mouse_position(2));
	 %    gsize = size(params.grid);

	 %    igrid = params.grid + params.padding;
	 %    xc = params.mouse_position(1);
	 %    yc = params.mouse_position(2);
	 %    [jj,ii] = ind2sub( [gsize(2), gsize(1)], ...
		% params.control_point_selected );
	 %    xa = min(jj+1, gsize(2));
	 %    xb = max(jj-1, 1);
		% ya = min(ii+1, gsize(1));
		% yb = max(ii-1, 1);
		% kxb = (ii-1)*gsize(2) + xb;
		% kyb = (yb-1)*gsize(2) + jj;

		% if jj == gsize(2) - 2
		% 	set(params.handle_grid_hlines(k), ...
	 %    	'XData', [], ...
	 %    	'YData', [] )
	 %    else
		%     set(params.handle_grid_hlines(k), ...
		%     	'XData', [xc igrid(ii, xa, 2)], ...
		%     	'YData', [yc igrid(ii, xa, 1)] )
		% end

		% if jj == 2
		% 	set(params.handle_grid_hlines(kxb), ...
	 %    	'XData', [], ...
	 %    	'YData', [] )
		% else
		%     set(params.handle_grid_hlines(kxb), ...
		%     	'XData', [igrid(ii, xb, 2) xc], ...
		%     	'YData', [igrid(ii, xb, 1) yc] )
		% end

		% if ii == gsize(1) - 2
		% 	set(params.handle_grid_vlines(k), ...
	 %    	'XData', [], ...
	 %    	'YData', [] )
	 %   	else
		%     set(params.handle_grid_vlines(k), ...
		%     	'XData', [xc igrid(ya, jj, 2)], ...
		%     	'YData', [yc igrid(ya, jj, 1)] )
		% end

		% if ii == 2
		% 	set(params.handle_grid_vlines(kyb), ...
	 %    	'XData', [], ...
	 %    	'YData', [] )
		% else
		%     set(params.handle_grid_vlines(kyb), ...
		%     	'XData', [igrid(yb, jj, 2) xc], ...
		%     	'YData', [igrid(yb, jj, 1) yc] )
		% end

		% imsize = size(params.imt);
		% x1 = params.padding + 1; 
		% x2 = params.padding + imsize(1);
		% y1 = params.padding + 1; 
		% y2 = params.padding + imsize(2);
		% params.imsgrid = get_checkerboard(igrid, params.imsgrid);
		% %params.imsgrid = params.imsgrid + .5*params.imch;
		% set(params.handle_imgrid, 'CData', params.imsgrid);
	end
	set_params(params);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
	% hObject    handle to figure1 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	disp('figure1_WindowButtonUpFcn');
	params = get_params(); 
	if( isempty(params) ), return, end
	if(params.empty), return, end
	if (strcmpi(params.option, 'apply_grid')), return, end

	if(ishandle(params.handle_control_point_selected))
	    set(params.handle_control_point_selected, ...
	    	'MarkerSize', 8,  'MarkerFaceColor', 'none');
	    update_grid();
	end
	
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
	% hObject    handle to figure1 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Hint: delete(hObject) closes the figure
	params = get_params();
	if( isempty(params) ), return, end

	params.optvideo.pause = true;
	set_params(params);
		
	% Stop the timer function
	if (params.empty)
		% Close the Window
		delete(hObject);
	else
		if isfield(params, 'gui_timer')
			stop(params.gui_timer); 
			delete(params.gui_timer);
		end
		% Close the Window
		delete(hObject);
	end
	
	close all;

% --------------------------------------------------------------
function saves_image_Callback(hObject, eventdata, handles)
	% hObject    handle to saves_image (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params(); 
	if(params.empty), return, end

	[filename, pathname] = uiputfile({'*.png';'*.jpg';'*.bmp'}, ...
		'Save the deformed image as');
	if isequal(filename,0)
		disp('User selected Cancel')
		return
	end
	imwrite(params.imt, fullfile(pathname, filename));

% --------------------------------------------------------------
function grid_Callback(hObject, eventdata, handles)
	% hObject    handle to grid (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------
function resize_grid_Callback(hObject, eventdata, handles)
	% hObject    handle to refine_grid (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if( params.empty ), return, end
	if (strcmpi(params.option, 'apply_grid')), return, end

	imsize = size(params.im);

	prompt = {...
		'Enter the number of control points in dimesion 1:', ...
		'Enter the number of control points in dimesion 2:' };
	dlgtitle = 'Grid size';
	dims = [1 35];
	definput = {...
		num2str(floor(imsize(1) / params.spacing(1) + 1)), ...
		num2str(floor(imsize(1) / params.spacing(2) + 1))};

	stop(params.gui_timer);
	answer = inputdlg(prompt,dlgtitle,dims,definput);
	start(params.gui_timer);

	if( isempty(answer) ), return, end	

	gpoints = zeros(2,1);
	gpoints(1) = str2double(answer{1}) - 1;
	gpoints(2) = str2double(answer{2}) - 1;

	params.spacing = zeros(2,1);
	params.spacing(1) = floor(imsize(1) / gpoints(1));
	params.spacing(2) = floor(imsize(2) / gpoints(2));
	params.uniform = make_init_grid( params.spacing, ...
		imsize(1:2) );
	params.grid = params.uniform;
	params.imt = params.im;

	if params.isrgb
		switch(params.channel)
			case 0
		        params.uniformrgb = params.uniform;
		    	params.gridrgb = params.grid;
		    case 1
		    	params.uniformred = params.uniform;
		    	params.gridred = params.grid;
		    case 2
		    	params.uniformgreen = params.uniform;
		    	params.gridgreen = params.grid;
		    case 3
		        params.uniformblue = params.uniform;
		    	params.gridblue = params.grid;
		    otherwise 
		end
	end

	set_params(params);

	figure(params.f)
	set_params(params);
	figure(params.handles.figure1)

	show_slices();
	show_grid();


% --------------------------------------------------------------
function save_grid_Callback(hObject, eventdata, handles)
	% hObject    handle to save_grid (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params(); 
	if( params.empty ), return, end
	if (strcmpi(params.option, 'apply_grid')), return, end

	[filename, pathname] = uiputfile({'*.mat'}, ...
		'Save the grid of control points as');

	imsize = size(params.im);

	gridcp = params.grid;
	spacing = params.spacing;

	gridcp(:,:,1) = gridcp(:,:,1) / imsize(1);
	gridcp(:,:,2) = gridcp(:,:,2) / imsize(2);

	
	spacing(1) = spacing(1) / imsize(1);
	spacing(2) = spacing(2) / imsize(2);

	save(fullfile(pathname, filename), ...
		'gridcp', 'spacing');


% --------------------------------------------------------------
function load_grid_Callback(hObject, eventdata, handles)
	% hObject    handle to load_grid (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if( params.empty ), return, end
	if ~params.optvideo.pause
		params.optvideo.pause = true;
		set_params(params);
	end
	
	params.optvideo.pause = false;

	[file, fpath] = uigetfile( ...
		{'*.mat'; '*.*'}, 'Load a grid' );
	if isequal(file,0)
		disp('User selected Cancel')
		return
	end
	
	load(fullfile(fpath, file));

	if params.optvideo.status
		imsize = [params.video.Height, params.video.Width];
	else
		imsize = size(params.im);
	end

	spacing(1) = floor(imsize(1) * spacing(1));
	spacing(2) = floor(imsize(2) * spacing(2));
	gridcp(:,:,1) = floor(imsize(1) * gridcp(:,:,1));
	gridcp(:,:,2) = floor(imsize(2) * gridcp(:,:,2));

	params.spacing = spacing;
	params.grid = gridcp;

	params.uniform = make_init_grid( params.spacing, ...
		[imsize(1:2)] );

	set_params(params);

	try
		figure(params.f);
		set_params(params);
		figure(params.handles.figure1)
	catch
		disp('Second figure dont exist');
	end

	if params.optvideo.status
		%params.video.CurrentTime = 0.0;
		show_video();
	else
		bspline_tranformix_image();
		if strcmpi(params.option, 'build_grid')
			show_slices();
			show_grid();
		elseif strcmpi(params.option, 'apply_grid')
			disp('apply_grid');
			params = get_params();
			set(params.handle_imgrid, 'CData', params.imt);
		end
	end


% --------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
	% hObject    handle to options (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------
function build_grid_Callback(hObject, eventdata, handles)
	% hObject    handle to build_grid (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	params.option = 'build_grid';
	set_params(params);
	set_menu_checks();

	if ( params.empty )
		return
	end
	
	clear_figure1();

	params.empty = true;
	params.optvideo.status = false;
	params.optvideo.pause = true;
	if isfield(params, 'gui_timer')
		params = rmfield(params, 'gui_timer');
	end
	set_params(params);


% --------------------------------------------------------------
function apply_grid_Callback(hObject, eventdata, handles)
	% hObject    handle to apply_grid (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	%if(isempty(params)), return, end
	params.option = 'apply_grid';
	set_params(params);
	set_menu_checks();
	if (params.empty)
		return
	end

	clear_figure1();

	params.empty = true;
	params.optvideo.status = false;
	params.optvideo.pause = true;
	if isfield(params, 'gui_timer')
		params = rmfield(params, 'gui_timer');
	end
	set_params(params);
	
	


% --------------------------------------------------------------
function open_video_Callback(hObject, eventdata, handles)
	% hObject    handle to open_video (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	disp('open_video');

	params = get_params();
	if ( strcmpi(params.option, 'build_grid') )
		warndlg(...
			'You must activate the "Apply a Grid" option', ...
			'Inavild Option');
		return
	end

	[file, fpath] = uigetfile( ...
		{'*.avi';'*.mp4';'*.*'}, ...
		'Open a video' );
	if isequal(file,0)
		disp('User selected Cancel')
		return
	end
	params.empty = false;
	params.optvideo.status = true;
	params.optvideo.pause = false;
	params.video = VideoReader(fullfile(fpath, file));
	vsize = [params.video.Height, params.video.Width];

	% Get the power of max refinement steps
	powermax = min( floor( log2(vsize*.25) ) );

	% Get b-spline grid spacing in x and y direction
	params.spacing = [2^powermax, 2^powermax];
	params.uniform = make_init_grid( params.spacing, ...
		vsize );
	params.grid = params.uniform;
	params.handle_grid_axes = params.handles.axes1;
	params.handle_control_point_selected = [];
	params.handle_control_points = [];

	params.frame = readFrame(params.video);
	params.framet = params.frame;
	params.handle_imgrid = imshow(params.frame, 'Parent', ...
		params.handle_grid_axes);
	%params.handle_imgrid = image(params.frame, 'Parent', ...
	%	params.handle_grid_axes);

	set_params(params);
	show_video()


% --------------------------------------------------------------
function video_Callback(hObject, eventdata, handles)
	% hObject    handle to video (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------
function reset_video_Callback(hObject, eventdata, handles)
	% hObject    handle to reset_video (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if ( ~params.optvideo.status )
		return
	end
	params.video.CurrentTime = 0.0;

	if params.optvideo.pause
		params.optvideo.pause = false;
		set_params(params);
		show_video();
	else
		set_params(params);
	end

% --------------------------------------------------------------
function play_video_Callback(hObject, eventdata, handles)
	% hObject    handle to play_video (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if ( ~params.optvideo.status )
		return
	end
	if params.optvideo.pause
		params.optvideo.pause = false;
		set_params(params);
		show_video();
	end
	


% --------------------------------------------------------------
function pause_video_Callback(hObject, eventdata, handles)
	% hObject    handle to pause_video (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if ( ~params.optvideo.status )
		return
	end
	params.optvideo.pause = true;
	set_params(params);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
	% hObject    handle to slider1 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Hints: get(hObject,'Value') returns position of slider
	%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
	disp('slider1_CreateFcn')
	params = get_params();
	params.alpha = get(hObject,'Value');
	if( isempty(params) ), return, end
	if(params.empty)
		set_params(params);		
		return
	end
	%disp(params.alpha)

	imsize = size(params.imt);
	x1 = params.padding + 1; 
	x2 = params.padding + imsize(1);
	y1 = params.padding + 1; 
	y2 = params.padding + imsize(2);

	params.imsgrid( x1:x2, y1:y2, : ) = params.imt;

	igrid = params.grid + params.padding;
	params.imsgrid = get_checkerboard(igrid, params.imsgrid, params.alpha);
	set(params.handle_imgrid, 'CData', params.imsgrid);	

	set_params(params);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
	% hObject    handle to slider1 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    empty - handles not created until after all CreateFcns called

	% Hint: slider controls usually have a light gray background.
	if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	    set(hObject,'BackgroundColor',[.9 .9 .9]);
	end


% --------------------------------------------------------------
function rgb_Callback(hObject, eventdata, handles)
	% hObject    handle to rgb (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------
function red_Callback(hObject, eventdata, handles)
	% hObject    handle to red (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if( isempty(params) ), return, end
	if( params.empty ), return, end
	if( ~params.isrgb ), return, end

	set(params.handles.red,'Checked','on');
	set(params.handles.green,'Checked','off');
	set(params.handles.blue,'Checked','off');
	set(params.handles.full,'Checked','off');

	params.channel = 1;

	if strcmpi(params.option, 'apply_grid')
		if ~params.optvideo.status
			params.im = params.imrgb;
			params.im(:,:,2:3) = 0.0;
			set_params(params);
			bspline_tranformix_image();

			params = get_params();
			set(params.handle_imgrid, 'CData', params.imt);
		end
		set_params(params);
	else
		params.im = params.imred;
		params.uniform = params.uniformred;
		params.grid = params.gridred;
		params.spacing = params.spacingred;

		figure(params.f);
		set_params(params);
		figure(params.handles.figure1)
		set_params(params);

		bspline_tranformix_image();
		show_slices();
		show_grid();
	end

% --------------------------------------------------------------
function green_Callback(hObject, eventdata, handles)
	% hObject    handle to green (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if( isempty(params) ), return, end
	if( params.empty ), return, end
	if( ~params.isrgb ), return, end

	set(params.handles.red,'Checked','off');
	set(params.handles.green,'Checked','on');
	set(params.handles.blue,'Checked','off');
	set(params.handles.full,'Checked','off');

	params.channel = 2;

	if strcmpi(params.option, 'apply_grid')
		if ~params.optvideo.status
			params.im = params.imrgb;
			params.im(:,:,1:2:3) = 0.0;
			set_params(params);
			bspline_tranformix_image();

			params = get_params();
			set(params.handle_imgrid, 'CData', params.imt);
		end
		set_params(params);
	else
		params.im = params.imgreen;
		params.uniform = params.uniformgreen;
		params.grid = params.gridgreen;
		params.spacing = params.spacinggreen;

		figure(params.f);
		set_params(params);
		figure(params.handles.figure1)
		set_params(params);

		bspline_tranformix_image();
		show_slices();
		show_grid();
	end

% --------------------------------------------------------------
function blue_Callback(hObject, eventdata, handles)
	% hObject    handle to blue (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if( isempty(params) ), return, end
	if( params.empty ), return, end
	if( ~params.isrgb ), return, end

	set(params.handles.red,'Checked','off');
	set(params.handles.green,'Checked','off');
	set(params.handles.blue,'Checked','on');
	set(params.handles.full,'Checked','off');

	params.channel = 3;

	if strcmpi(params.option, 'apply_grid')
		if ~params.optvideo.status
			params.im = params.imrgb;
			params.im(:,:,1:2) = 0.0;
			set_params(params);
			bspline_tranformix_image();

			params = get_params();
			set(params.handle_imgrid, 'CData', params.imt);
		end
		set_params(params);
	else
		params.im = params.imblue;
		params.uniform = params.uniformblue;
		params.grid = params.gridblue;
		params.spacing = params.spacingblue;

		set_params(params);
		figure(params.f);
		set_params(params);
		figure(params.handles.figure1)

		bspline_tranformix_image();
		show_slices();
		show_grid();
	end

% --------------------------------------------------------------
function full_Callback(hObject, eventdata, handles)
	% hObject    handle to full (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	params = get_params();
	if( isempty(params) ), return, end
	if( params.empty ), return, end
	if( ~params.isrgb ), return, end

	set(params.handles.red,'Checked','off');
	set(params.handles.blue,'Checked','off');
	set(params.handles.green,'Checked','off');
	set(params.handles.full,'Checked','on');

	params.channel = 0;

	if strcmpi(params.option, 'apply_grid')
		if ~params.optvideo.status
			params.im = params.imrgb;
			set_params(params);
			bspline_tranformix_image();

			params = get_params();
			set(params.handle_imgrid, 'CData', params.imt);
		end
		set_params(params);
	else
		params.im = params.imrgb;
		params.uniform = params.uniformrgb;
		params.grid = params.gridrgb;
		params.spacing = params.spacingrgb;

		figure(params.f)
		set_params(params);
		figure(params.handles.figure1)
		set_params(params);

		bspline_tranformix_image();
		show_slices();
		show_grid();
	end