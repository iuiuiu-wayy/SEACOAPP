from flask import Flask
from dash import Dash
from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from leftmenu import Left_menu, ids, menu_logics
from Figure import call_figure, rightMenu


# class Config:
#     #database
#     SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:project@localhost/ProjectDB'
#     SQLALCHEMY_ECHO = False
#     SQLALCHEMY_TRACK_MODIFICATIONS = False

app = Flask(__name__, instance_relative_config=True)
# app.config.from_object(Config)

dash_app = Dash(
    server=app, 
    suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.SANDSTONE],
    url_base_pathname='/')

dash_app.layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col(html.H1('SEA-COAPP'), 
            align='center', width=3)
        ], justify='center'),
        dbc.Row([
            dbc.Col(Left_menu, width=3,), #style={'background-color': 'coral'}),
            dbc.Col(rightMenu, width=9,), #style={'background-color': 'red'})
            # dbc.Col(tabs, width=True)
        ])     
    ])
    
])
dash_app.ids = ids


menu_logics(dash_app)


call_figure(dash_app)


if __name__ == '__main__':
    app.run(debug=True)