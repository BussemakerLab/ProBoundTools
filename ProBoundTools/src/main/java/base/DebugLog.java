package base;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

//import main.SELEX;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

public class DebugLog
{
    public static boolean DEBUG = true;
    public final static SimpleDateFormat format2=new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss");
    public final static SimpleDateFormat format=new SimpleDateFormat("HH:mm:ss.SSS");
    private static PrintWriter defaultLogger;
    private static String defaultLogFolder ;
    
    public static void log(Object message)
    {
    	String msg=formatLogMessage(message);
    	if(defaultLogFolder!=null)
    	{
        	getDefaultLogger().println(msg);
        	getDefaultLogger().flush();
    	}
        if (DEBUG)
        {
        	System.out.println(msg);
        }
    }
    
    public static void setDefaultLogFolder(String folder)
    {
    	defaultLogFolder =  folder;
    }
    
    public static String formatLogMessage(Object message)
    {
    	if(message instanceof Throwable)
    	{
    		Throwable e=(Throwable)message;
    		e.printStackTrace();
    	}
        String fullClassName = Thread.currentThread().getStackTrace()[3].getClassName();            
        String className = fullClassName.substring(fullClassName.lastIndexOf(".") + 1);
        String methodName = Thread.currentThread().getStackTrace()[3].getMethodName();
        int lineNumber = Thread.currentThread().getStackTrace()[3].getLineNumber();
        
       return format.format(new Date())+":"+ className + "." + methodName + "():" + lineNumber+ "#" + Thread.currentThread().getName() +"\t"+message;
    
    }
    
    private static PrintWriter getDefaultLogger()
    {
    	if(defaultLogger!=null)
    		return defaultLogger;
    	try
		{
    		String path=defaultLogFolder+"/"+format2.format(new Date())+".log";
			defaultLogger = new PrintWriter(path);
		} catch (FileNotFoundException e)
		{
			throw new RuntimeException(e);
		}
    	return defaultLogger;
    }
    
    public static void verbose(String trueFalse)
    {
    	DEBUG = Boolean.parseBoolean(trueFalse);
    }
    
    public static void verbose(boolean trueFalse)
    {
    	DEBUG =trueFalse;
    }
}