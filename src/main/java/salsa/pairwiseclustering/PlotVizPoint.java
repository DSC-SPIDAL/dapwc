package salsa.pairwiseclustering;

public class PlotVizPoint
{
	private double _x;
	private double _y;
	private double _z;
	private int _index;
	private int _cluster;
	private String _label;

	public PlotVizPoint(double x, double y, double z, int index, int cluster)
	{
		_x = x;
		_y = y;
		_z = z;
		_index = index;
		_cluster = cluster;
	}

	public PlotVizPoint(double x, double y, double z, int index, int cluster, String label)
	{
		_x = x;
		_y = y;
		_z = z;
		_index = index;
		_cluster = cluster;
		_label = label;
	}

	public final double getX()
	{
		return _x;
	}
	public final void setX(double value)
	{
		_x = value;
	}

	public final double getY()
	{
		return _y;
	}
	public final void setY(double value)
	{
		_y = value;
	}

	public final double getZ()
	{
		return _z;
	}
	public final void setZ(double value)
	{
		_z = value;
	}

	public final int getIndex()
	{
		return _index;
	}
	public final void setIndex(int value)
	{
		_index = value;
	}

	public final int getCluster()
	{
		return _cluster;
	}
	public final void setCluster(int value)
	{
		_cluster = value;
	}

	public final String getLabel()
	{
		return _label;
	}
	public final void setLabel(String value)
	{
		_label = value;
	}

    // TODO - fix plotviz point
	/*public final XElement ToPvizPointElement()
	{
	   return new XElement("point", new XElement("key", _index), new XElement("clusterkey", _cluster), new XElement("label", _label), new XElement("location", new XAttribute("x", _x), new XAttribute("y", _y), new XAttribute("z", _z)));
	}*/

	public final PlotVizPoint Clone()
	{
		return new PlotVizPoint(_x, _y, _z, _index, _cluster, _label);
	}
}