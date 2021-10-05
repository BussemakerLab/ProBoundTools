package base;

import java.util.Comparator;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

public class SequenceComparator implements Comparator<CountObject>
{
	//@Override
	public int compare(CountObject arg0, CountObject arg1)
	{
		return arg0.compareToByKey(arg1);
	}
}
