from django.conf.urls.defaults import patterns, include, url
from django.conf import settings

# Uncomment the next two lines to enable the admin:
#from django.contrib import admin
#admin.autodiscover()

urlpatterns = patterns('gui.views',
    # Examples:
    #url(r'^$', 'stanalyzer.views.home', name='home'),
    # url(r'^stanalyzer/', include('stanalyzer.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^$', 'index'),
    url(r'^help/', 'help'),
    url(r'^login/', 'login'),
    url(r'^logout/', 'logout'),
    url(r'^desktop/', 'desktop'),
    url(r'^project/', 'prjView'),
    url(r'^project_new/', 'prjView_new'),
    url(r'^project_update/', 'prjView_update'),
    url(r'^project_delete/', 'prjView_delete'),
    url(r'^user/', 'userView'),
    url(r'^stanalyzer/', 'stanalyzer'),
    url(r'^get_flist/', 'stanalyzer_info'),
    url(r'^get_structure/', 'stanalyzer_info'),
    url(r'^get_segment/', 'stanalyzer_info'),
    url(r'^getdb/', 'getDBinfo'),
    url(r'^path_validation/', 'pathValidation'),
    url(r'^permission_write/', 'info_permission_write'),
    url(r'^permission_exec/', 'info_permission_exec'),
    url(r'^sendJob/', 'stanalyzer_sendJob'),
    url(r'^job_view/', 'jobView'),
    url(r'^jobView_jqGrid_prj/', 'jobView_jqGrid_prj'),
    url(r'^jobView_jqGrid_job/', 'jobView_jqGrid_job'),
    url(r'^jobView_jqGrid_para/', 'jobView_jqGrid_para'),
    url(r'^usrView_jqGrid_usr/', 'usrView_jqGrid_usr'),
    url(r'^resultView_jqGrid_results/', 'resultView_jqGrid_results'),
    url(r'^resultView_del_results/', 'resultView_jqGrid_del_results'),
    url(r'^resultView_jqGrid/', 'resultView_jqGrid'),
    url(r'^create_user/', 'usrView_jqGrid_create_user'),
    url(r'^update_user/', 'usrView_jqGrid_create_user'),
    url(r'^del_user/', 'usrView_jqGrid_del_user'),
    url(r'^download/', 'makeDownload'),
    url(r'^resultView_DBmanager/', 'resultView_DBmanager'),
    
    url(r'^mediaLink/(?P<file_name>.*$)', 'mediaLink'),
    url(r'^wys_FileManager/', 'wysFileManager'),
    
    url(r'^toydb_data/', 'toyView_data'),
#    url(r'^toydb_prj/', 'toyView_prj'),
#    url(r'^toydb_usr/', 'toyView_usr'),
#    url(r'^set/root', 'set_root'),
#    url(r'^session/(?P<page_name>\w{0,50})/$', 'chksession'),
)


