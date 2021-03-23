import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State


Left_menu = html.Div([
    ##### 4 lat and lon options
    dbc.Row([
        dbc.Col(html.Label(['Lat S:']),width=3),
        dbc.Col(dbc.Input(type="number", min=-90, max=90, id='lats', value=-7.125, step=0.001),width=9)
    ]),
    dbc.Row([
        dbc.Col(html.Label(['Lat N:']),width=3),
        dbc.Col(dbc.Input(type="number", min=-90, max=90, id='latn', value=-5.875, step=0.001),width=9)
    ]),
    dbc.Row([
        dbc.Col(html.Label(['Lon L:']),width=3),
        dbc.Col(dbc.Input(type="number", min=0, max=360, id='lonl', value=108.875, step=0.001),width=9)
    ]),
    dbc.Row([
        dbc.Col(html.Label(['Lon R:']),width=3),
        dbc.Col(dbc.Input(type="number", min=0, max=360, id='lonr', value=110.125, step=0.001),width=9)
    ]),    

    ##### option model only or model + correction 
    dbc.Row([
        dbc.Col(html.Label(['Mode ']),width=3),
        dbc.Col(
            dcc.Dropdown(
                    id='mode',
                    options=[
                        {'label': 'Model', 'value':'model'},
                        {'label': 'Model + correction', 'value':'modelcorrection'},
                    ],
                    value = 'model',
                    multi = False,
                    clearable=False,
                    # style = {'width':"50%"},
                ),width=9)
    ]),

    ##### radio item for montly or yearly
    dbc.Row([
        dbc.Col(html.Label(['Option']),width=3),
        dbc.Col(
            dbc.RadioItems(
            options=[
                {"label": "Monthly", "value": 'monthly'},
                {"label": "Yearly", "value": 'yearly'},
            ],
            value='yearly',
            id="options",
            inline=True,
        ),width=9)
    ]),

    ###### pilih tahun koreksi  TODO havent choose the range correctly
    dbc.Row([
        dbc.Col(html.Label(['Select the year of correction:']),width=12),
    ]),
    dbc.Row([
        dbc.Col(html.Label(['YStart:']), width=3),
        dbc.Col(dbc.Input(type="number", min=2006, max=2018, id='YStart', value=2006, step=1),width=9)
    ]),       
    dbc.Row([
        dbc.Col(html.Label(['YEnd:']), width=3),
        dbc.Col(dbc.Input(type="number", min=2006, max=2018, id='YEnd', value=2018, step=1),width=9)
    ]),


    ###### pilih tahun koreksi  TODO havent choose the range correctly
    dbc.Row([
        dbc.Col(html.Label(['Select the year of output:']),width=12),
    ]),
    dbc.Row([
        dbc.Col(html.Label(['YStart:']), width=3),
        dbc.Col(dbc.Input(type="number", min=2006, max=2040, id='YOutStart', value=2006, step=1),width=9)
    ]),       
    dbc.Row([
        dbc.Col(html.Label(['YEnd:']), width=3),
        dbc.Col(dbc.Input(type="number", min=2006, max=2040, id='YOutEnd', value=2040, step=1),width=9)
    ]),

    ###### pilih variabel
    dbc.Row([
        dbc.Col(html.Label(['Pilih variabel: ']),width=3),
        dbc.Col(
            dcc.Dropdown(
                    id='var',
                    options=[
                        {'label': 'SSH', 'value':'SSH'},
                        {'label': 'SST', 'value':'SST'},
                        {'label': 'Salt', 'value':'SALT'},
                        {'label': 'U Velocity', 'value':'U_VEL'},
                        {'label': 'V Velocity', 'value':'V_VEL'},
                    ],
                    value = 'SSH',
                    multi = False,
                    clearable=False,
                    # style = {'width':"50%"},
                ),width=9)
    ]),

    ##### the two button
    dbc.Row([
        dbc.Col(dbc.Button("DOWNLOAD DATA", outline=True, color="primary", className="mr-1", id='button_plot', style = {'width':"100%"}), width=6),
        dbc.Col(dbc.Button("EXPORT", outline=True, color="primary", className="mr-1", id='button_export',style = {'width':"100%"}), width=6),
    ])
])

ids = ['lats', 'latn', 'lonl', 'lonr', 'mode', "options", 'YStart', 'YEnd', 'YOutStart', 'YOutEnd', 'var']

def menu_logics(app):
    
    @app.callback(
        [Output('mode', 'value'),
        Output('mode', 'options')],
        [Input('var', 'value')]
    )
    def mode_logic(var):
        options=[
            {'label': 'Model', 'value':'model'},
            {'label': 'Model + correction', 'value':'modelcorrection'},
        ]
        value = 'modelcorrection'
        if not var == 'SSH':
            options = [
                {'label': 'Model', 'value':'model'}
            ]
            value = 'model'
        return value, options