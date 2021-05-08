import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash_extensions import Download
from dash_extensions.snippets import send_file
from netCDF4 import date2num
import netCDF4
import numpy as np
import pandas as pd
from scipy import interpolate

ExportWindow = html.Div([
    dbc.Modal(
        [
            dbc.ModalHeader("Header"),
            dbc.ModalBody(
                "Change the backdrop of this modal with the radio buttons"
            ),
            dbc.ModalFooter(
                dbc.Button(
                    "Close", id="close_window", className="ml-auto"
                )
            ),
        ],
        id="export_Window",
        centered=True,
        backdrop='static',
    ),
])

def Export_window_Callbacks(app):
    @app.callback(
        Output("export_Window", "is_open"),
        [Input("button_export", "n_clicks"), Input("close_window", "n_clicks")],
        [State("export_Window", "is_open")],
    )
    def toggle_export(n1, n2, is_open):
        if n1 or n2:
            return not is_open
        return is_open

    @app.callback(
        Output('export_Window', 'children'),
        [Input('button_export', 'n_clicks')]
    )
    def dynamic_content(n1):
        if not app.downloadable:    
            content= [
                dbc.ModalHeader("No data avaiable"),
                dbc.ModalBody(
                    "Please download the data first before exporting it"
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close", id="close_window", className="ml-auto"
                    )
                ),
            ]
        else:
            content = [
                dbc.ModalHeader("Select data to download"),
                dbc.ModalBody([
                    Download(id="download_ti"),
                    Download(id='download_clim'),
                    Download(id='download_aave'),
                    Download(id='download_anom'),
                    Download(id='download_anom_csv'),
                    Download(id='download_data'),
                    dbc.Row([
                        dbc.Col(width=3),
                        dbc.Col([
                            dbc.ButtonGroup(
                                [
                                    dbc.Button("Trend and intercept (nc)", id='button_ti'),
                                    dbc.Button("Climatology (nc)", id='button_clim'),
                                    dbc.Button("Areal average of data (csv)", id='button_aave'),
                                    dbc.Button("Anomaly (nc)", id='button_anom'),
                                    dbc.Button("Anomaly (csv)", id='button_anom_csv'),
                                    dbc.Button("Data (nc)", id='button_data'),                            
                                ],
                                vertical=True,
                            )
                        ])
                        
                    ])
                    
                ]
                    # "Please download the data first before exporting it"
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close", id="close_window", className="ml-auto"
                    )
                ),
            ]
            
        return content
    
    ######### callback for trend and intercepts
    @app.callback(
        Output('download_ti', 'data'),
        [Input('button_ti', 'n_clicks')],
        prevent_initial_call=True
    )
    def click_dl_ti(n_clicks):
        ncfile = netCDF4.Dataset('output_ti.nc',mode='w',format='NETCDF4_CLASSIC')
        lat_dim = ncfile.createDimension('lat', len(app.latitudes))     
        lon_dim = ncfile.createDimension('lon', len(app.longitudes))

        ncfile.title='My model data'

        lat = ncfile.createVariable('lat', np.float64, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float64, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'

        var1 = ncfile.createVariable('Intercept',np.float64,('lat','lon')) # note: unlimited dimension is leftmost
        var1.units = '-' # degrees Kelvin
        var1.standard_name = 'Intercept' # this is a CF standard name
        var1.long_name = 'Intercept of model'

        var2 = ncfile.createVariable('Trend',np.float64,('lat','lon')) # note: unlimited dimension is leftmost
        var2.units = '-' # degrees Kelvin
        var2.standard_name = 'Trend' # this is a CF standard name
        var2.long_name = 'Trend of model'

        nlats = len(lat_dim)
        nlons = len(lon_dim)

        lat[:]=app.latitudes
        lon[:]=app.longitudes
        print('res shape is ', app.res_b.shape)
        print('lon lan', nlons, nlats)
        var1[:,:] = app.res_b
        var2[:,:] = app.res_a
        ncfile.close()
        return send_file(
            'output_ti.nc', filename='output ti.nc'
        )


    ######## callback for climatology
    @app.callback(
        Output('download_clim', 'data'),
        [Input('button_clim', 'n_clicks')],
        prevent_initial_call=True
    )
    def click_dl_clim(n_clicks):
        ncfile = netCDF4.Dataset('output_clim.nc',mode='w',format='NETCDF4_CLASSIC')
        lat_dim = ncfile.createDimension('lat', len(app.latitudes))     
        lon_dim = ncfile.createDimension('lon', len(app.longitudes))

        ncfile.title='My model data'

        lat = ncfile.createVariable('lat', np.float64, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float64, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'

        if app.options == 'monthly':
            month_dim = ncfile.createDimension('month', 12)
            month = ncfile.createVariable('month', np.int8, ('month',))
            month.units = '0-11'
            month.long_name = 'index of month'

            var1 = ncfile.createVariable('Climatology',np.float64,('month','lat','lon')) # note: unlimited dimension is leftmost
            var1.units = 'Monthly climatology' # degrees Kelvin
            var1.standard_name = 'Climatology' # this is a CF standard name
            var1.long_name = 'Climatogology of model'
            
            lat[:]=app.latitudes
            lon[:]=app.longitudes
            month[:]=[i for i in range(12)]
            print('monthlyclim', app.spatial_monthly_clim.shape)
            print(12,len(lat), len(lon))
            var1[:,:,:] = app.spatial_monthly_clim

        else:
            var1 = ncfile.createVariable('Climatology',np.float64,('lat','lon')) # note: unlimited dimension is leftmost
            var1.units = 'Yearly climatology' # degrees Kelvin
            var1.standard_name = 'Climatology' # this is a CF standard name
            var1.long_name = 'Climatogology of model'

            lat[:]=app.latitudes
            lon[:]=app.longitudes
            var1[:,:] = app.agg_data.mean(0)
        nlats = len(lat_dim)
        nlons = len(lon_dim)


        ncfile.close()
        return send_file(
            'output_clim.nc', filename='output clim.nc'
        )

    ######## callback for aave csv
    @app.callback(
        Output('download_aave', 'data'),
        [Input('button_aave', 'n_clicks')],
        prevent_initial_call=True
    )
    def click_dl_aave(n_clicks):
        outputfilename = 'output_aave.csv'
        outdf = pd.DataFrame({
            'Date': app.agg_label,
            'Data': app.agg_data.mean(1).mean(1)
        })
        outdf.to_csv(outputfilename,index=False)
        return send_file(
            'output_aave.csv', filename='output_aave.csv'
        )

    ######## callback for anomaly nc
    @app.callback(
        Output('download_anom', 'data'),
        [Input('button_anom', 'n_clicks')],
        prevent_initial_call=True
    )
    def click_dl_anom(n_clicks):
        ncfile = netCDF4.Dataset('output_anomaly.nc',mode='w',format='NETCDF4_CLASSIC')
        lat_dim = ncfile.createDimension('lat', len(app.latitudes))     
        lon_dim = ncfile.createDimension('lon', len(app.longitudes))
        time_dim = ncfile.createDimension('time', None)


        lat = ncfile.createVariable('lat', np.float64, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float64, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        
        time_data = ncfile.createVariable('time', np.float64, ('time',) )
        time_data.long_name = 'time'

        anom = ncfile.createVariable('anom',np.float64,('time','lat','lon'))
        anom.units = app.var_units
        
        ncfile.title='My model data'

        dummy_latitudes = app.obs_latitudes
        dummy_longitudes = app.obs_longitudes
        def z_func(z):
            return z
        if len(app.obs_longitudes)==1:
            dummy_longitudes = app.obs_longitudes.tolist() + [app.obs_longitudes[0] + 0.001]
            def z_func(z):
                temp_z = np.append(z, z, axis=1)
                return temp_z
        if len(app.obs_latitudes) == 1:
            dummy_latitudes = app.obs_latitudes.tolist() + [app.obs_latitudes[0] + 0.001]
            def z_func(z):
                temp_z = np.append(z, z, axis=0)
                return temp_z
            
        if app.options == 'yearly':
            time_data.units = 'hours since 1800-01-01'
            anom.long_name = 'Yearly anomaly'
            lat[:]=app.latitudes
            lon[:]=app.longitudes
            tmp_clim = app.agg_obs_data.mean(0)
            if app.mode == 'modelcorrection':


                f = interpolate.interp2d(dummy_longitudes, dummy_latitudes, z_func(tmp_clim), kind='linear')
                
                clim_regrid = f(app.longitudes, app.latitudes)
                spatial_anomaly = app.agg_data - clim_regrid
                
            else:
                spatial_anomaly = app.agg_data - app.agg_data.mean(0)
            times = date2num(app.agg_label, time_data.units)



        else:
            time_data.units = 'hours since 1800-01-01'
            anom.long_name = 'Monthly anomaly'
            lat[:]=app.latitudes
            lon[:]=app.longitudes
            spatial_anomaly = []

            if app.mode == 'modelcorrection':
                regrid_spatial_montly_clim = []
                for i in range(12):
                    # print('z shape', z_func(app.base_spatial_montly_clim[i]).shape)
                    # print('len lon', len(dummy_latitudes), len(dummy_longitudes))
                    # print('app.obs_latitudes ', dummy_latitudes)
                    # print('app.obs_longitudes ', dummy_longitudes)
                    f = interpolate.interp2d(dummy_longitudes, dummy_latitudes,  z_func(app.base_spatial_montly_clim[i]), kind='linear')
                    regrid_spatial_montly_clim.append(f(app.longitudes, app.latitudes))
                
                for i in range(app.agg_data.shape[0]):
                    tmp_map = app.agg_data[i,:,:] - regrid_spatial_montly_clim[i%12]
                    spatial_anomaly.append(tmp_map)
        
            else:
                for i in range(app.agg_data.shape[0]):
                    tmp_map = app.agg_data[i,:,:] - app.spatial_monthly_clim[i%12]
                    spatial_anomaly.append(tmp_map)
            
            spatial_anomaly = np.array(spatial_anomaly)
            times = date2num(app.agg_label, time_data.units)
        
        anom[:,:,:] = spatial_anomaly
        
        time_data[:]=times

        # var1 = ncfile.createVariable('anomaly',np.float64,('time', 'lat', 'lon')) # note: unlimited dimension is leftmost
        # var1.units = '-' # degrees Kelvin
        # var1.standard_name = 'Intercept' # this is a CF standard name
        # var1.long_name = 'Intercept of model'

        # var2 = ncfile.createVariable('Trend',np.float64,('lat','lon')) # note: unlimited dimension is leftmost
        # var2.units = '-' # degrees Kelvin
        # var2.standard_name = 'Trend' # this is a CF standard name
        # var2.long_name = 'Trend of model'

        # nlats = len(lat_dim)
        # nlons = len(lon_dim)

        # lat[:]=app.latitudes
        # lon[:]=app.longitudes
        # print('res shape is ', app.res_b.shape)
        # print('lon lan', nlons, nlats)
        # var1[:,:] = app.res_b
        # var2[:,:] = app.res_a
        ncfile.close()
        return send_file(
            'output_anomaly.nc', filename='output_anomaly.nc'
        )

    ######## callback for anomaly csv
    @app.callback(
        Output('download_anom_csv', 'data'),
        [Input('button_anom_csv', 'n_clicks')],
        prevent_initial_call=True
    )
    def click_dl_anom_csv(n_clicks):
        dummy_latitudes = app.obs_latitudes
        dummy_longitudes = app.obs_longitudes
        def z_func(z):
            return z
        if len(app.obs_longitudes)==1:
            dummy_longitudes = app.obs_longitudes.tolist() + [app.obs_longitudes[0] + 0.001]
            def z_func(z):
                temp_z = np.append(z, z, axis=1)
                return temp_z
        if len(app.obs_latitudes) == 1:
            dummy_latitudes = app.obs_latitudes.tolist() + [app.obs_latitudes[0] + 0.001]
            def z_func(z):
                temp_z = np.append(z, z, axis=0)
                return temp_z

        if app.options == 'yearly':
            tmp_clim = app.agg_obs_data.mean(0)
            if app.mode == 'modelcorrection':
                f = interpolate.interp2d(dummy_longitudes, dummy_latitudes, z_func(tmp_clim), kind='linear')    
                clim_regrid = f(app.longitudes, app.latitudes)
                spatial_anomaly = app.agg_data - clim_regrid
            else:
                spatial_anomaly = app.agg_data - app.agg_data.mean(0)
            
            

        else:
            spatial_anomaly = []
            if app.mode == 'modelcorrection':
                regrid_spatial_montly_clim = []
                for i in range(12):
                    f = interpolate.interp2d(dummy_longitudes, dummy_latitudes,  z_func(app.base_spatial_montly_clim[i]), kind='linear')
                    regrid_spatial_montly_clim.append(f(app.longitudes, app.latitudes))
                
                for i in range(app.agg_data.shape[0]):
                    tmp_map = app.agg_data[i,:,:] - regrid_spatial_montly_clim[i%12]
                    spatial_anomaly.append(tmp_map)
        
            else:
                for i in range(app.agg_data.shape[0]):
                    tmp_map = app.agg_data[i,:,:] - app.spatial_monthly_clim[i%12]
                    spatial_anomaly.append(tmp_map)
            
            spatial_anomaly = np.array(spatial_anomaly)
        
        outdf = pd.DataFrame({
            'Date': app.agg_label,
            'Data': spatial_anomaly.mean(1).mean(1)
        })
        outputfilename = 'output_anomaly.csv'
        outdf.to_csv(outputfilename,index=False)
        return send_file(
            'output_anomaly.csv', filename='output_anomaly.csv'
        )


    ######## callback for data nc
    @app.callback(
        Output('download_data', 'data'),
        [Input('button_data', 'n_clicks')],
        prevent_initial_call=True
    )
    def click_dl_data(nclicks):

        ncfile = netCDF4.Dataset('output_data.nc',mode='w',format='NETCDF4_CLASSIC')
        lat_dim = ncfile.createDimension('lat', len(app.latitudes))     
        lon_dim = ncfile.createDimension('lon', len(app.longitudes))
        time_dim = ncfile.createDimension('time', None)


        lat = ncfile.createVariable('lat', np.float64, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float64, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        
        time_data = ncfile.createVariable('time', np.float64, ('time',) )
        time_data.long_name = 'time'

        anom = ncfile.createVariable('data',np.float64,('time','lat','lon'))
        anom.units = app.var_units
        
        ncfile.title='My model data'

        time_data.units = 'hours since 1800-01-01'
        anom.long_name = 'Spatial data'
        lat[:]=app.latitudes
        lon[:]=app.longitudes

        times = date2num(app.agg_label, time_data.units)

        anom[:,:,:] = app.agg_data
        
        time_data[:]=times

        app.agg_data

        ncfile.close()
        return send_file(
            'output_data.nc', filename='output_data.nc'
        )