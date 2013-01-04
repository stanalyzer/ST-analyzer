from gui.models import User, Project, Job, Parameter
from django.contrib import admin

admin.site.register(User);

# Defining project form    
class ProjectAdmin (admin.ModelAdmin):
    fieldsets = [
        (None, {'fields':['user', 'name']}),
        ('Date Information', {'fields':['date']}),
        ('Default PBS script', {'fields':['pbs']}),
        ]
    list_filter = ['date'] 
    date_hierarchy = 'date'
        
admin.site.register(Project, ProjectAdmin);

# Defining Job form
class JobAdmin (admin.ModelAdmin):
    fieldsets = [
        (None, {'fields':['proj', 'name']}),
        ('Analyzer', {'fields':['anaz']}),
        ('Date Information', {'fields':['stime', 'etime']}),
        ('Results', {'fields':['status', 'output']}),
        ]
    list_filter = ['stime']
    date_hierarchy = 'etime'

admin.site.register(Job, JobAdmin);

admin.site.register(Parameter);
#admin.site.register(Analyzer);
