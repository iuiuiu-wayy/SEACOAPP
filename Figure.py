import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from datetime import datetime
import numpy as np
import netCDF4
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from dateutil.rrule import rrule, MONTHLY, YEARLY
from calendar import monthrange


rightMenu = dcc.Loading(
    type="circle", 
    fullscreen=False,
    children = html.Div([
        html.Div(id='div_id'),
        dcc.Graph(id='figure_id')
    ])
)    


colorscale= px.colors.sequential.Plasma
# colorscale[0] = "#FFFFFF00"

def daily2morthlyyearly(options, start,end,data):
    startID = 0
    resultData = []
    resultLable = []
    if options == 'yearly':
        iterator = YEARLY
    else:
        iterator = MONTHLY
    for m in rrule(iterator,dtstart=datetime(start, 1,1), until=datetime(end,12,31)):
        if not options == 'yearly':
            n_days = monthrange(m.year,m.month)[1]
        else:
            if m.year % 4 == 0:
                n_days = 366
            else:
                n_days = 365
        tmp = data[startID:startID+n_days, :, :].mean(0)
        resultData.append(tmp)
        resultLable.append(m)
        startID = startID + n_days
    
    return resultLable, resultData


def call_figure(app):
    @app.callback(
        [Output('figure_id', 'figure'),
        Output('div_id', 'children')],
        [Input('button_plot', 'n_clicks')],
        [State(idss, 'value') for idss in app.ids]
    )

    def create_figure(n_clicks, *args):
        if n_clicks is None:
            raise PreventUpdate

        app.downloadable = False
        # print(args)
        ######## modelll
        CLIMCHANGE_INPUT = "http://tides.big.go.id:8080/thredds/dodsC/InaROMS/InaROMS_rcp45_2006_2040.nc"
        dataset = netCDF4.Dataset(CLIMCHANGE_INPUT)
        latitudes = dataset.variables['LAT'][:] 
        longitudes = dataset.variables['LON'][:]
        lats, latn, lonl, lonr, mode, options, YStart, YEnd, YOutStart, YOutEnd, var = args
        app.var_units = dataset[var].units
        app.var_long_name= dataset[var].long_name
        app.args = args
        app.mode = mode
        app.options = options

        regionLon = [ lonl , lonr ]
        regionLat = [lats, latn]
        time = dataset.variables['OCEAN_TIME'][:]
        time_conver = netCDF4.num2date(time,dataset.variables['OCEAN_TIME'].units, only_use_cftime_datetimes=False)

        time_li = time_conver.tolist().index((datetime(int(YStart),1,1,12)))
        time_ui = time_conver.tolist().index((datetime(int(YEnd),12,31, 12)))

        latli = np.argmin( np.abs( latitudes - regionLat[0] ) )
        latui = np.argmin( np.abs( latitudes - regionLat[1] ) )

        # longitude lower and upper index
        lonli = np.argmin( np.abs( longitudes - regionLon[0] ) )
        lonui = np.argmin( np.abs( longitudes - regionLon[1] ) )  

        if var == 'SSH':
            data = dataset.variables[var][ : , latli:latui , lonli:lonui ]
        else:
            data = dataset.variables[var][ : , : , latli:latui , lonli:lonui ]
            data = data[:,0,:,:]
        
        if mode == 'modelcorrection':
            model_fac = data[time_li:time_ui].mean(0).mean(0).mean(0)

            #####obs for SSH
            obs_file = 'sealevel_indo.nc'
            obs_ds = netCDF4.Dataset(obs_file)

            time_obs = obs_ds['time'][:]
            time_converted = netCDF4.num2date(time_obs,obs_ds.variables['time'].units, only_use_cftime_datetimes=False)
            lats_obs = obs_ds['latitude'][:]
            lons_obs = obs_ds['longitude'][:]


            time_li_obs = time_converted.tolist().index((datetime(int(YOutStart),1,1)))
            time_ui_obs = time_converted.tolist().index((datetime(int(YEnd),12,31)))

            # latitude lower and upper index
            latli_obs = np.argmin( np.abs( lats_obs - regionLat[0] ) )
            latui_obs = np.argmin( np.abs( lats_obs - regionLat[1] ) )

            # longitude lower and upper index
            lonli_obs = np.argmin( np.abs( lons_obs - regionLon[0] ) )
            lonui_obs = np.argmin( np.abs( lons_obs - regionLon[1] ) )
            obsdata = obs_ds['adt'][time_li:time_ui, latli_obs:latui_obs, lonli_obs:lonui_obs]
            obs_fac = obsdata.mean(0).mean(0).mean(0)


            ##### calculate correction factor
            correction_factor = obs_fac - model_fac
            # print('correction factor is ', correction_factor)

            # generate data for model
            data = data + correction_factor

            agg_obs_label, agg_obs_data = daily2morthlyyearly(options, int(YStart), int(YEnd), obsdata)
            app.agg_obs_data = np.array(agg_obs_data)
            obs_lat, obs_lon = np.meshgrid(lons_obs[lonli_obs:lonui_obs],lats_obs[latli_obs:latui_obs])

            app.obs_latitudes = lats_obs[latli_obs:latui_obs]
            app.obs_longitudes = lons_obs[lonli_obs:lonui_obs]

            app.agg_obs_label = agg_obs_label
            app.obs_lat = obs_lat.flatten().data
            app.obs_lon = obs_lon.flatten().data
            app.obs_cent_lon = lons_obs[lonli_obs:lonui_obs].mean()
            app.obs_cent_lat = lats_obs[latli_obs:latui_obs].mean()
        

        #generate the desired output
        agg_label, agg_data = daily2morthlyyearly(options, int(YOutStart), int(YOutEnd), data)
        agg_data = np.array(agg_data)
        app.agg_data = agg_data  



        lon, lat = np.meshgrid(longitudes[lonli:lonui],latitudes[latli:latui])
        app.longitudes=longitudes[lonli:lonui]
        app.latitudes=latitudes[latli:latui]

        lon = lon.flatten().data
        lat = lat.flatten().data

        app.agg_label = agg_label
        app.agg_data = agg_data
        app.lon = lon
        app.lat = lat

        app.cent_lat = latitudes[latli:latui].mean()
        app.cent_lon = longitudes[lonli:lonui].mean()
        # plot_array = agg_data[0,:,:].flatten()


        y = app.agg_data.mean(1).mean(1)
        x = app.agg_label

        fig = go.Figure(data=go.Scatter(
            x=x,
            y=y,
            mode='lines+markers',
        ))
        # x = ["%04i-%02i" % (m.year, m.month) for m in agg_label]
        # print(plot_array)
        # fig = go.Figure(
        #     go.Scattermapbox(
        #         lat=lat,
        #         lon = lon,
        #         mode='markers',
        #         text=plot_array,
        #         marker=go.scattermapbox.Marker(
        #             cmax = plot_array.max(),
        #             cmin = plot_array.min(),
        #             color = plot_array,
        #             colorscale = colorscale
        #             # size=9
        #         ),
        #     )
        # )

        # fig.update_layout(
        #     margin=dict(t=0,b=0,r=0,l=0),
        #     autosize=True,
        #     hovermode='closest',
        #     mapbox_style="carto-positron",
        #     mapbox=dict(
        #         # accesstoken=mapbox_access_token,
        #         bearing=0,
        #         center=dict(
        #             lat=app.cent_lat,
        #             lon=app.cent_lon
        #         ),
        #         pitch=0,
        #         zoom=5
        #     ),
        # )
        


        # datafig = [
        #     go.Scattermapbox(
        #         lon=lon,
        #         lat=lat,
        #         mode='markers',
        #         text=plot_array,
        #         marker=go.scattermapbox.Marker(
        #             # cmax=2.5,
        #             # cmin=-2.5,
        #             color=plot_array,
        #             # colorscale=colorscale
        #         ),
        #     )
        # ]

        # layout = go.Layout(
        #     margin=dict(t=0,b=0,r=0,l=0),
        #     autosize=True,
        #     hovermode='closest',
        #     showlegend=False,
        #     mapbox=dict(
        #         bearing=0,
        #         center=dict(
        #             lat=latitudes.mean(),
        #             lon=longitudes.mean()
        #         ),
        #         pitch=0,
        #         zoom=0,
        #         style='dark'
        #     ),
        # )
        # fig = dict(data=data, layout=layout)

        if mode == 'modelcorrection':

            div = dcc.Dropdown(
                id='right_dropdown',
                options=[
                    {'label':'Trend', 'value':'trend'},
                    {'label':'Intercept', 'value':'intercept'},
                    {'label':'First data', 'value':'data'},
                    {'label':'Base climatology', 'value':'base_clim'},
                    {'label':'Model climatology', 'value':'clim'},
                    {'label': 'Areal average scatter plot', 'value':'scatter'},
                    {'label': 'Anomaly', 'value':'anomaly'}
                ],
                style = {'width':"50%"},
                value= 'scatter',
                multi = False,
                clearable=False,
            )
        else:
            div = dcc.Dropdown(
                id='right_dropdown',
                options=[
                    {'label':'Trend', 'value':'trend'},
                    {'label':'Intercept', 'value':'intercept'},
                    {'label':'First data', 'value':'data'},
                    {'label':'Mocel climatology', 'value':'clim'},
                    {'label': 'Areal average scatter plot', 'value':'scatter'},
                    {'label': 'Anomaly', 'value':'anomaly'}
                ],
                style = {'width':"50%"},
                value= 'scatter',
                multi = False,
                clearable=False,
            )
        app.second_call = False
        app.downloadable=True
        return fig, div


    @app.callback(
        Output('figure_id', 'figure'),
        [Input('right_dropdown', 'value')]
    )
    def second_fig_call(right_dropdown):
        print('entering second calll')
        if not app.second_call:
            n_t, n_i, n_j = app.agg_data.shape
            res_a = []
            res_b = []

            for i in range(n_i):
                tmp_a = []
                tmp_b = []

                for j in range(n_j):
                    y = app.agg_data[:,i,j]
                    x = np.array(range(len(y)))
                    a, b = np.polyfit(x, y, 1)
                    tmp_a.append(a)
                    tmp_b.append(b)
                res_a.append(tmp_a)
                res_b.append(tmp_b)

            app.res_a = np.array(res_a)
            app.res_b = np.array(res_b)

            # clim = app.agg_data.mean(0)
        if app.options == 'monthly':
            ######## calculate monthly climatology
            raw_montly = app.agg_data.mean(1).mean(1)
            monthly_clim = list()
            spatial_monthly_clim = list()
            
            for i in range(12):
                tmp_data = []
                for j in range(i, len(raw_montly), 12):
                    
                    tmp_data.append(raw_montly[j])
                tmp_data = np.array(tmp_data)
                monthly_clim.append(tmp_data.mean())

                spatial_tmp = list()
                # print('agg_data', app.agg_data.shape)
                # print('i index data', app.agg_data.shape[0])
                # test_tmp = list()
                n_index = app.agg_data.shape[0]
                # makii = app.agg_data[0,:,:]
                # test_tmp.append(makii)
                # test_tmp.append(makii)
                # print('test_tmp', test_tmp)
                # print('test_tmp_array', np.array(test_tmp).shape)
                # print('sample data luar',makii)
                for j2 in range(i, n_index, 12):
                    # print('setelah masuk loop')
                    # print('j ini', j2)
                # makii = [[0,0],[0,0]]
                    makii = app.agg_data[j,:,:]
                    spatial_tmp.append(makii)
                # spatial_tmp.append()
                # print(n_index, 'n index')
                # print('sample data',makii)
                spatial_tmp = np.array(spatial_tmp)
                
                spatial_monthly_clim.append(spatial_tmp.mean(0))
            #     print('spatial_tmp', spatial_tmp)
            # print('spatial_montly_clim', len(spatial_monthly_clim), spatial_monthly_clim[-1].shape)
            # print('agg_data shpae', app.agg_data.shape)

            monthly_clim = np.array(monthly_clim)
            app.monthly_clim = monthly_clim
            app.spatial_monthly_clim = np.array(spatial_monthly_clim)


            if app.mode == 'modelcorrection':
                raw_base_montly = app.agg_obs_data.mean(1).mean(1)
                base_monthly_clim = []
                base_spatial_montly_clim = list()
                for i in range(12):
                    base_tmp_data = []
                    base_tmp_spatial = []
                    for j in range(i,len(raw_base_montly), 12):
                        base_tmp_data.append(raw_base_montly[j])
                        base_tmp_spatial.append(app.agg_obs_data[j,:,:])
                    base_tmp_data = np.array(base_tmp_data)
                    base_monthly_clim.append(base_tmp_data.mean())

                    base_tmp_spatial = np.array(base_tmp_spatial)
                    base_spatial_montly_clim.append(base_tmp_spatial.mean(0))
                base_monthly_clim = np.array(base_monthly_clim)
                app.base_spatial_montly_clim = np.array(base_spatial_montly_clim)
                # print('random check ', app.base_spatial_montly_clim.shape)

        if right_dropdown == 'trend':
            plot_array = app.res_a.flatten()
            lat_used = app.lat
            lon_used = app.lon
            app.lineplot = False
        elif right_dropdown == 'intercept':
            plot_array = app.res_b.flatten()
            lat_used = app.lat
            lon_used = app.lon
            app.lineplot = False
        elif right_dropdown == 'clim':

            if app.options == 'yearly':
                plot_array = app.agg_data.mean(0).flatten()
                lat_used = app.lat
                lon_used = app.lon
                app.lineplot = False
            else:
                app.lineplot = True
                x = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct','Nov', 'Dec' ]
                y = monthly_clim
        elif right_dropdown == 'base_clim':
            if app.options == 'yearly':
                plot_array = app.agg_obs_data.mean(0).flatten()
                lat_used = app.obs_lat
                lon_used = app.obs_lon
                app.lineplot = False
            else:
                app.lineplot = True
                x = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct','Nov', 'Dec' ]
                y = base_monthly_clim
                    
        elif right_dropdown == 'scatter':
            app.lineplot = True
            y = app.agg_data.mean(1).mean(1)
            x = app.agg_label
        elif right_dropdown == 'anomaly':
            app.lineplot = True
            if app.mode == 'modelcorrection':
                x = app.agg_obs_label
                if app.options == 'yearly':
                    tmp_clim = app.agg_obs_data.mean(0).mean(0).mean(0)
                    y = app.agg_data.mean(1).mean(1) - tmp_clim
                else:
                    anomaly_montly = []
                    for i in range(n_t):
                        anomaly_montly.append(raw_montly[i] - base_monthly_clim[i%12])
                    y = np.array(anomaly_montly)
            else:                    
                x = app.agg_label
                if app.options == 'yearly':
                    tmp_clim = app.agg_data.mean(0).mean(0).mean(0)
                    y = app.agg_data.mean(1).mean(1) - tmp_clim
                else:
                    anomaly_montly =  []
                    for i in range(n_t):
                        anomaly_montly.append(raw_montly[i] - monthly_clim[i%12])
                    y = np.array(anomaly_montly)
        else:
            plot_array = app.agg_data[2,:,:].flatten()
            lat_used = app.lat
            lon_used = app.lon
            app.lineplot = False

        if app.lineplot:
            
            fig = go.Figure(data=go.Scatter(
                x=x,
                y=y,
                mode='lines+markers',
            ))
        else:

            #### re arraange data mask
            new_lat = []
            new_lon = []
            new_data = []

            for i in range(len(plot_array)):
                if not plot_array[i] == 0:
                    new_data.append(plot_array[i])
                    new_lon.append(lon_used[i])
                    new_lat.append(lat_used[i])
            new_data = np.array(new_data)
            new_lon = np.array(new_lon)
            new_lat = np.array(new_lat)
            # print(plot_array)
            fig = go.Figure(
                go.Scattermapbox(
                    lat=new_lat, #app.lat,
                    lon = new_lon, #app.lon,
                    mode='markers',
                    text= new_data, #plot_array,
                    marker=go.scattermapbox.Marker(
                        cmax = new_data.max(), #plot_array.max(),
                        cmin = new_data.min(), #plot_array.min(),
                        color = new_data,
                        colorscale = colorscale,
                        size=13
                    ),
                )
            )

            fig.update_layout(
                margin=dict(t=0,b=0,r=0,l=0),
                autosize=True,
                hovermode='closest',
                mapbox_style="carto-positron",
                mapbox=dict(
                    # accesstoken=mapbox_access_token,
                    bearing=0,
                    center=dict(
                        lat=app.cent_lat,
                        lon=app.cent_lon
                    ),
                    pitch=0,
                    zoom=5
                ),
            )
        app.downloadable=True
        return fig
