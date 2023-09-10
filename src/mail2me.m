function mail2me(subject,content,path)
MailAddress = 'a1714170867@163.com';
password = 'mikutjw419027';  
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','smtp.163.com');
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
sendmail('a1714170867@163.com',subject,content,path);