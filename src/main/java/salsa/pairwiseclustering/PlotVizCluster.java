package salsa.pairwiseclustering;

public class PlotVizCluster
{
	private int _clusterNumber;
	private String _label;

	public PlotVizCluster(int clusterNumber, String label)
	{
		_clusterNumber = clusterNumber;
		_label = label;
	}

	public PlotVizCluster(int clusterNumber)
	{
		_clusterNumber = clusterNumber;
	}

	public PlotVizCluster(String label)
	{
		_label = label;
	}

	public final int getClusterNumber()
	{
		return _clusterNumber;
	}
	public final void setClusterNumber(int value)
	{
		_clusterNumber = value;
	}

	public final String getLabel()
	{
		return _label;
	}
	public final void setLabel(String value)
	{
		_label = value;
	}

    // TODO - fix plotviz cluster
	/*public final XElement ToClusterElement(Color color, boolean isDefault, double size)
	{
		return new XElement("cluster", new XElement("key", _clusterNumber), new XElement("label", _label), new XElement("visible", 1), new XElement("default", isDefault ? 1 : 0), new XElement("color", new XAttribute("r", color.R), new XAttribute("g", color.G), new XAttribute("b", color.B), new XAttribute("a", color.A)), new XElement("size", size));
	}*/
}