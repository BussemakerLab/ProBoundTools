package base;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.lang.reflect.Method;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Properties;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

public class Util
{
	public static int atoi(String base, int start)
	{
		int result = 0;
		for (int i = start; i < base.length(); i++)
		{
			char c = base.charAt(i);
			int v = c - '0';
			if (v >= 0 && v < 10)
			{
				result = result * 10 + v;
			} else
			{
				break;
			}
		}
		return result;
	}

	private static final double LOG2 = Math.log(2);

	public static double log2(double a)
	{
		return Math.log(a) / LOG2;
	}

	public static void printMemeoryUsage()
	{
		int mb = 1024 * 1024;
		// Getting the runtime reference from system
		Runtime runtime = Runtime.getRuntime();
		// Print used memory
		DebugLog.log("Used Memory: "
				+ (runtime.totalMemory() - runtime.freeMemory()) / mb + " MB");
		// Print free memory
		DebugLog.log("Free Memory: " + runtime.freeMemory() / mb + " MB");
		DebugLog.log("Max Memory:  " + runtime.maxMemory() / mb + " MB");
	}
	
	/**
	 * Called in R.
	 * @return
	 */
	public static Object[] getMemeoryUsage()
	{
		int mb = 1024 * 1024;

		Runtime runtime = Runtime.getRuntime();
		String[] attributes=new String[]{"Free Memory","Used Memory","Max Memory"};
		String[] values=new String[]{
				( (runtime.totalMemory() - runtime.freeMemory()) /mb)+" MB",
				( runtime.freeMemory()/mb)+" MB",
				( runtime.maxMemory()/mb)+" MB"};
		
		return new Object[]{attributes,values };
	}

	public static File createTempDirectory() throws Exception
	{
		File temp;

		temp = File.createTempFile("temp", Long.toString(System.nanoTime()));

		if (!(temp.delete()))
		{
			throw new IOException("Could not delete temp file: "
					+ temp.getAbsolutePath());
		}

		if (!(temp.mkdir()))
		{
			throw new IOException("Could not create temp directory: "
					+ temp.getAbsolutePath());
		}

		return (temp);
	}
/*
	public static SequencingRunInfo cloneWithoutSample(SequencingRunInfo oldInfo)
	{
		try
		{
			SequencingRunInfo SequencingRunInfo = (SequencingRunInfo) cloneObjectBySetterGetter(oldInfo);
			return SequencingRunInfo;
		}
		catch(Exception ex)
		{
			DebugLog.log(ex);
			throw new RuntimeException(ex);
		}
	}

	public static Sample cloneSample(Sample oldInfo)
	{
		try
		{
			Sample sample = (Sample) cloneObjectBySetterGetter(oldInfo);
			return sample;
		}
		catch(Exception ex)
		{
			DebugLog.log(ex);
			throw new RuntimeException(ex);
		}
	}*/

	@SuppressWarnings("unchecked")
	public static Object cloneObjectBySetterGetter(Object objectIn) throws Exception
	{
		Class cls = objectIn.getClass();
		Object clonedObj = objectIn.getClass().newInstance();
		Method[] gettersAndSetters = cls.getMethods();

		for (int i = 0; i < gettersAndSetters.length; i++)
		{
			String methodName = gettersAndSetters[i].getName();
			try
			{
				if (methodName.startsWith("get"))
				{
					cls.getMethod(methodName.replaceFirst("get", "set"),
									gettersAndSetters[i].getReturnType())
							.invoke(clonedObj,
									gettersAndSetters[i].invoke(objectIn));
				} else if (methodName.startsWith("is"))
				{
					cls.getMethod(methodName.replaceFirst("is", "set"),
									gettersAndSetters[i].getReturnType())
							.invoke(clonedObj,
									gettersAndSetters[i].invoke(objectIn));
				}

			} catch (NoSuchMethodException e)
			{
			} catch (IllegalArgumentException e)
			{
			}

		}

		return clonedObj;
	}
	

	public static Properties loadProperties(String configFile)
	{
		DebugLog.log("Reading configuration file:" + configFile);
		Properties prop = new Properties();
		FileInputStream input;
		try
		{
			input = new FileInputStream(configFile);
			prop.load(input);
		} catch (Exception e)
		{
			DebugLog.log(e);
			DebugLog.log(e); throw new RuntimeException("Can't load the config file:"
					+ configFile);
		}
		return prop;
	}
	
	public static boolean isAbsolutePath(String path)
	{
		if(path.startsWith("/")) //unix
		{
			return true;
		}
		else if(path.length()>2 && path.charAt(1)==':') //windows C:\\
		{
			return true;
		}
		return false;
	}
	
	public static String getMD5(String text)
	{
        try {
            MessageDigest md = MessageDigest.getInstance("MD5");
            //Using MessageDigest update() method to provide input
            byte[] buffer = text.getBytes();
            int numOfBytesRead = buffer.length;
            md.update(buffer, 0, numOfBytesRead);
            byte[] hash = md.digest();
            return bytesToHex(hash);
            
        } catch (Exception ex) {
            DebugLog.log(ex);
            throw new RuntimeException(ex);
        }
          
    }

	//copy from http://stackoverflow.com/questions/9655181/convert-from-byte-array-to-hex-string-in-java
	final protected static char[] hexArray = "0123456789ABCDEF".toCharArray();
	public static String bytesToHex(byte[] bytes) {
	    char[] hexChars = new char[bytes.length * 2];
	    for ( int j = 0; j < bytes.length; j++ ) {
	        int v = bytes[j] & 0xFF;
	        hexChars[j * 2] = hexArray[v >>> 4];
	        hexChars[j * 2 + 1] = hexArray[v & 0x0F];
	    }
	    return new String(hexChars);
	}
	
	public static boolean isEmpty(String s)
	{
		return s==null || s.trim().length()==0;
	}


	public static void radixSort(CountObject[] list, int size)
	{
		int len = list[0].getKey().getLength();
		int BUCKET_SIZE = 6;
		LinkedList<CountObject>[] buckets = new LinkedList[ (int)Math.pow(4, BUCKET_SIZE) ];
		for(int j=0;j<buckets.length;j++)
		{
			buckets[j] = new LinkedList<CountObject>();
		}
		
		long mark = (long) (1 << (2*BUCKET_SIZE)) - 1;
		int shift = 0;

		int steps = (int) Math.ceil(len*1.0/BUCKET_SIZE);
		for(int i=0;i< steps;i++)
		{

			//System.out.println(Arrays.toString(list));
			
			for(int k=0;k<size;k++)
			{
				int idx = (int) ( ( list[k].getKey().getValue() & mark ) >> shift ) ;
				buckets[idx].add(list[k]);
			}
			
			int index = 0;
			for(int j=0;j<buckets.length;j++)
			{
				for(CountObject o:buckets[j])
				{
					list[index++]= o ;
				}
				buckets[j].clear();
			}
			
			mark <<= (BUCKET_SIZE*2);
			shift += (BUCKET_SIZE*2);

		}
	}
	
	public static int getNumberOfCores()
	{
		return Runtime.getRuntime().availableProcessors();
	}
	
	public static void outputColumnBasedArrays(PrintWriter writer, Object[] data, String[] header)
	{
		for(String h:header)
		{
			writer.print(h+"\t");
		}
		writer.println();
		
		int colNumber = data.length;
		int rowNumber = Array.getLength(data[0]);
		//DebugLog.log("Rows : "+ rowNumber) ; 
		
		for(int i=0;i<rowNumber; i++)
		{
			for(int j=0;j<colNumber ; j++)
			{
				writer.print( Array.get(data[j],i) + "\t");
			}
			writer.println();
		}
	}
	
	public static void output1DArray(PrintWriter writer, Object[] data)
	{
		for(Object h:data)
		{
			writer.print(h+"\t");
		}
		writer.println();
	}
	
	public static void output1DArray(PrintWriter writer, Collection data, int skip)
	{
		int i=0;
		for(Object h:data)
		{
			i++;
			if(i>skip)
			{
				writer.print(h+"\t");
			}
		}
		writer.println();
	}

	
	public static String joinStrings(Collection<String> strs, String seperator)
	{
		StringBuffer sb=new StringBuffer();
		for(String s:strs)
		{
			sb.append(s);
			sb.append(seperator);
		}
		sb.setLength(sb.length()-seperator.length());
		return sb.toString();
	}

	private static int MAX_NUM_OF_THREAD ;
	
	static
	{
		int cores = Util.getNumberOfCores();
		MAX_NUM_OF_THREAD = Math.min(cores, 4)*2;
	}
	
	public static int getMaxThreadNumber()
	{
		return MAX_NUM_OF_THREAD;
	}
	
	/**
	 * Called in R.
	 * @param num
	 */
	public static void setMaxThreadNumber(Integer num)
	{
		if(num>0)
		{
			MAX_NUM_OF_THREAD = num;
			DebugLog.log("Max number of thread : "+ num);
		}
	}
	
}
