from django.db import models

# Create your models here.
# default User name: admin, Password = 12345

class User (models.Model):
    uid  = models.CharField(primary_key=True, max_length=20);
    pwd = models.CharField(max_length=100);
    email = models.EmailField(null=True, blank=True);
    level = models.IntegerField();

class Project (models.Model):
    # Using primary_key with auto increment
    user = models.ForeignKey(User, on_delete=models.CASCADE);
    name = models.CharField(max_length=200);
    date  = models.DateTimeField('project created');
    pbs   = models.TextField();

#class Analyzer (models.Model):
    # Using primary_key with auto increment
#    name = models.CharField(max_length=200);
#    path = models.FilePathField();
class path_input (models.Model):
    # Using primary_key with auto increment
    proj = models.ForeignKey(Project, on_delete=models.CASCADE);
    path = models.CharField(max_length=300);
    
class path_output (models.Model):
    # Using primary_key with auto increment
    proj = models.ForeignKey(Project, on_delete=models.CASCADE);
    path = models.CharField(max_length=300);
    
class path_python (models.Model):
    # Using primary_key with auto increment
    proj = models.ForeignKey(Project, on_delete=models.CASCADE);
    path = models.CharField(max_length=300);
    
    
class Job (models.Model):
    # Using primary_key with auto increment
    name = models.CharField(max_length=200);
    proj = models.ForeignKey(Project, on_delete=models.CASCADE);
    #anaz = models.ForeignKey(Analyzer);
    anaz = models.TextField(null=True, blank=True);
    status = models.CharField(max_length=20);
    output = models.CharField(max_length=300);
    stime = models.DateTimeField();
    etime = models.DateTimeField(null=True, blank=True);

class Parameter (models.Model):
    # Using primary key with auto increment
    job = models.ForeignKey(Job, on_delete=models.CASCADE);
    anaz = models.CharField(max_length=100);
    para = models.CharField(max_length=100);
    val = models.CharField(max_length=100);
    status = models.CharField(max_length=20)

class Outputs (models.Model):
    # Using primary key with auto increment
    job   = models.ForeignKey(Job, on_delete=models.CASCADE);
    name  = models.CharField(max_length=100);
    img   = models.TextField(null=True, blank=True);
    txt   = models.TextField(null=True, blank=True);
    gzip = models.CharField(max_length=100);
