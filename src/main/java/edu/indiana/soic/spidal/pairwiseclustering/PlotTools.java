package edu.indiana.soic.spidal.pairwiseclustering;

public class PlotTools
{
    // TODO - fix PlotTools
/*	public static void CreatePlotWithCenters(String centerFile, String pointsFile, String clusterNumberFile, int numberOfCenterPointsToIncludeInEachCenterType, String centerPlotFile, String plotDescription)
	{
		*//* Generate all types of center clusters per cluster
		 * 
		 * Center clusters are,
		 *  1. Original Min Mean
		 *  2. MDS Min Mean
		 *  3. MDS Center of Gravity (CoG)
		 *  4. Overall Best
		 *  5. Bucket Fraction 0
		 *     Bucket Fraction 1 and so on
		 *     
		 * Number of center points to include in each center type = n 
		 * n <= N, which is the number of center points found for each center type by PWC
		 * N is specified through NumberOfCenters parameter in PWC
		 * 
		 * Assumes a center file from a PWC center finding run
		 * Assumes a points file, which has each point mapped to its cluster in the format 
		 *  PointNumber<TAB>Xcoord<TAB>Ycoord<TAB>Zcoord<TAB>ClusterNumber
		 *//*


		*//* Colors to use with PlotViz
		   reads color info from Matlab50.txt file *//*
		java.util.ArrayList<Color> matlab50Colors = GenerateMatlab50Colors();

		*//* XML elements to hold points and clusters to be used in PlotViz file *//*
		XElement clustersElement = new XElement("clusters");
		XElement pointsElement = new XElement("points");

		*//* Hashtable mapping point number to a PlotVizPoint data structure for the points in the given points file *//*
		java.util.Hashtable existingPointsTable = new java.util.Hashtable();

		*//* Maximum number of points int the points file *//*
		int maxpnum = 0;
		*//* Maximum number of clusters that points are mapped to in the points file*//*
		int maxcnum = 0;

		edu.indiana.soic.spidal.Boxspidal.general.Box<Integer> boxmaxpnum = new edu.indiana.soic.spidal.Boxspidal.general.Box<Integer>(maxpnum);
		edu.indiana.soic.spidal.generaloic.spidal.Box<Integer> boxmaxcnum = new edu.indiana.soic.spidal.Boxspidal.general.Box<Integer>(maxcnum);
		ProcessPointsFile(pointsFile, clusterNumberFile, clustersElement, pointsElement, boxmaxpnum, boxmaxcnum, existingPointsTable, matlab50Colors);
		maxpnum = boxmaxpnum.content;
		maxcnum = boxmaxcnum.content;

		*//* Table mapping each cluster (i.e. group) number to another table called method table
		 * method table maps each method (e.g. smallest distance mean, smallest MDS distance mean, etc.) name to the list center points for that particular method
		 * the order of points in the list is as same as in the given center file *//*
		java.util.Hashtable groupTable = ProcessCenterFile(centerFile);

		CreatePlotWithCentersInternal(centerPlotFile, plotDescription, clustersElement, pointsElement, maxpnum, existingPointsTable, maxcnum, matlab50Colors, groupTable, numberOfCenterPointsToIncludeInEachCenterType);

	}

	private static void CreatePlotWithCentersInternal(String centerPlotFile, String plotDescription, XElement clustersElement, XElement pointsElement, int maxpnum, java.util.Hashtable existingPointsTable, int maxcnum, java.util.ArrayList<Color> matlab50Colors, java.util.Hashtable groupTable, int numberOfCenterPointsToIncludeInEachCenterType)
	{
		++maxcnum;
		for (DictionaryEntry groupToMethodTable : groupTable)
		{
			int group = (int)groupToMethodTable.Key; // group is the original cluster number
			java.util.Hashtable methodTable = (java.util.Hashtable)groupToMethodTable.Value;
			int methodCount = methodTable.size();
			int tempCount = methodCount;
			for (DictionaryEntry methodToCenterPoints : methodTable)
			{
				String method = (String)methodToCenterPoints.Key; // method is one of smallest distance mean, smallest MDS mean, etc.

				// cluster number to be used in PlotViz for this center type 
				int methodNumber = methodCount - tempCount--;
				int clusterNumberForCenterType = group * methodCount + methodNumber + maxcnum;

				// cluster name to be used in PlotViz for this center type
				String centerTypeName = group + "" + method + ".centerpoints";

				// add an XML element to represent this center type as a cluster in PlotViz
				clustersElement.Add(CreateClusterElement(clusterNumberForCenterType, centerTypeName, matlab50Colors.get(group % matlab50Colors.size()), false, 2.0, methodNumber));

				java.util.ArrayList<CenterInfo> cps = (java.util.ArrayList<CenterInfo>)methodToCenterPoints.Value;
				// Picking the topmost n point for each method
				for (int i = 0; i < numberOfCenterPointsToIncludeInEachCenterType; i++)
				{
					CenterInfo cp = cps.get(i);
					PlotVizPoint p = (PlotVizPoint)existingPointsTable.get(cp.getPnum());
					pointsElement.Add(CreatePointElement(++maxpnum, clusterNumberForCenterType, ("cluster:" + group + "-idx:" + p.getIndex() + "method:" + method), p.getX(), p.getY(), p.getZ()));
				}
			}
		}

		XElement plotElement = CreatePlotElement(plotDescription, true);
		XElement plotvizElement = new XElement("plotviz");
		plotvizElement.Add(plotElement);
		plotvizElement.Add(clustersElement);
		plotvizElement.Add(pointsElement);
		plotvizElement.Save(centerPlotFile);

	}



	private static void ProcessPointsFile(String pointsFile, String clusterNumberFile, XElement clusters, XElement points, edu.indiana.soic.spidal.generaloic.spidal.Box<Integer> maxpnum, edu.indiana.soic.spidal.generaloic.spidal.Box<Integer> maxcnum, java.util.Hashtable pointsTable, java.util.ArrayList<Color> matlab50Colors)
	{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//		using (StreamReader preader = new StreamReader(pointsFile), creader = new StreamReader(clusterNumberFile))
		StreamReader preader = new StreamReader(pointsFile);
		StreamReader creader = new StreamReader(clusterNumberFile);
		try
		{
			java.util.HashSet<Integer> clusterNumbers = new java.util.HashSet<Integer>();
			maxpnum.content = -1;
			while (!preader.EndOfStream)
			{
				String pline = preader.ReadLine();
				String cline = creader.ReadLine();
				if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(pline) && !tangible.DotNetToJavaStringHelper.isNullOrEmpty(cline))
				{
					PlotVizPoint p = ReadPointLine(pline.trim());
					if (maxpnum.content < p.getIndex())
					{
						maxpnum.content = p.getIndex();
					}
					pointsTable.put(p.getIndex(), p);

					int cnum = ReadCnum(cline);
					p.setCluster(cnum);
					if (!clusterNumbers.contains(p.getCluster()))
					{
						clusterNumbers.add(p.getCluster());
						clusters.Add(CreateClusterElement(p.getCluster(), (new Integer(p.getCluster())).toString(CultureInfo.InvariantCulture), matlab50Colors.get(p.getCluster() % matlab50Colors.size()), true, 0.1, Glyphs.Hexagon2D));
					}
					points.Add(CreatePointElement(p.getIndex(), p.getCluster(), "", p.getX(), p.getY(), p.getZ()));
				}
			}
			maxcnum.content = clusterNumbers.Max();
		}
		finally
		{
			preader.dispose();
			creader.dispose();
		}
	}

	private static int ReadCnum(String line)
	{
		char[] sep = new char[] {' ', '\t'};
		String[] splits = line.split(sep, StringSplitOptions.RemoveEmptyEntries);
		return splits.length == 2 ? Integer.parseInt(splits[1]) : splits.length == 5 ? Integer.parseInt(splits[4]) : 0;
	}

	private static java.util.ArrayList<Color> GenerateMatlab50Colors()
	{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//		using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Salsa.PairwiseClusteringTPL.Matlab50.txt"))
		Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Salsa.PairwiseClusteringTPL.Matlab50.txt");
		try
		{
			if (stream != null)
			{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//				using (StreamReader reader = new StreamReader(stream))
				StreamReader reader = new StreamReader(stream);
				try
				{
					java.util.ArrayList<Color> colors = new java.util.ArrayList<Color>();
					char[] sep = new char[] {' ', '\t'};
					String[] splits;
					String split;
					int startIdx = 3;
					int r, g, b, a;
					while (!reader.EndOfStream)
					{
						String line = reader.ReadLine();
						if (!tangible.DotNetToJavaStringHelper.isNullOrEmpty(line))
						{
							splits = line.trim().split(java.util.regex.Pattern.quote(sep.toString()), -1);

							split = splits[0];
							r = Integer.parseInt(split.substring(startIdx, startIdx + (split.length() - (startIdx + 1))));

							split = splits[1];
							g = Integer.parseInt(split.substring(startIdx, startIdx + (split.length() - (startIdx + 1))));

							split = splits[2];
							b = Integer.parseInt(split.substring(startIdx, startIdx + (split.length() - (startIdx + 1))));

							split = splits[3];
							a = Integer.parseInt(split.substring(startIdx, startIdx + (split.length() - (startIdx + 1))));

							colors.add(Color.FromArgb(a, r, g, b));
						}
					}
					return colors;
				}
				finally
				{
					reader.dispose();
				}
			}
			else
			{
				throw new RuntimeException("Unable to load embedded resource: Matlab50.txt");
			}
		}
		finally
		{
			stream.dispose();
		}
	}

	private static PlotVizPoint ReadPointLine(String line)
	{
		char[] sep = new char[] {' ', '\t'};
		String[] splits = line.split(sep, StringSplitOptions.RemoveEmptyEntries);
		PlotVizPoint p = new PlotVizPoint(Double.parseDouble(splits[1]), Double.parseDouble(splits[2]), Double.parseDouble(splits[3]), Integer.parseInt(splits[0]), Integer.parseInt(splits[4]));
		return p;
	}

	private static CenterInfo ReadCenterLine(String line)
	{
		char[] sep = new char[] {' ', '\t'};
		char[] eqsep = new char[] {'='};
		String[] splits = line.split(sep, StringSplitOptions.RemoveEmptyEntries);
		int pnum = Integer.parseInt(splits[0].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1]);
		double measure = Double.parseDouble(splits[1].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1]);
		int methodIdx = 2;
		String source = "";
		double count = 0.0;
		if (splits[2].startsWith("Count"))
		{
			methodIdx = 4;
			count = Double.parseDouble(splits[2].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1]);
			source = splits[3].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1];
		}
		String method = splits[methodIdx].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1];
		int group = Integer.parseInt(splits[methodIdx + 1].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1]);
		String seqName = splits[methodIdx + 2].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1];
		for (int i = methodIdx + 3; i < splits.length - 4; ++i)
		{
			seqName += (" " + splits[i]);
		}
		int seqLength = Integer.parseInt(splits[splits.length - 4].split(java.util.regex.Pattern.quote(eqsep.toString()), -1)[1]);
		return new CenterInfo(pnum, measure, method, group, seqName, seqLength, source, count);
	}

	private static java.util.Hashtable ProcessCenterFile(String centerFile)
	{
//C# TO JAVA CONVERTER NOTE: The following 'using' block is replaced by its Java equivalent:
//		using (StreamReader reader = new StreamReader(centerFile))
		StreamReader reader = new StreamReader(centerFile);
		try
		{
			java.util.Hashtable groupTable = new java.util.Hashtable();
			while (!reader.EndOfStream)
			{
				CenterInfo cp = ReadCenterLine(reader.ReadLine());
				AddToGroupTable(groupTable, cp);
			}
			return groupTable;
		}
		finally
		{
			reader.dispose();
		}
	}

	private static void AddToGroupTable(java.util.Hashtable groupTable, CenterInfo cp)
	{
		if (groupTable.containsKey(cp.getCluster()))
		{
			java.util.Hashtable methodTable = (java.util.Hashtable)groupTable.get(cp.getCluster());
			if (methodTable.containsKey(cp.getMethod()))
			{
				// Need a list to maintain the order of points
				java.util.ArrayList<CenterInfo> cps = (java.util.ArrayList<CenterInfo>)methodTable.get(cp.getMethod());
				cps.add(cp);
			}
			else
			{
				// Need a list to maintain the order of points
				java.util.ArrayList<CenterInfo> cps = new java.util.ArrayList<CenterInfo>(java.util.Arrays.asList(new CenterInfo[] {cp}));
				methodTable.put(cp.getMethod(), cps);
			}
		}
		else
		{
			// Need a list to maintain the order of points
			java.util.ArrayList<CenterInfo> cps = new java.util.ArrayList<CenterInfo>(java.util.Arrays.asList(new CenterInfo[] {cp}));
			java.util.Hashtable methodTable = new java.util.Hashtable();
			methodTable.put(cp.getMethod(), cps);
			groupTable.put(cp.getCluster(), methodTable);
		}
	}


	private static XElement CreatePlotElement(String name, boolean glyphVisible)
	{
		XElement plot = new XElement("plot", new XElement("title", name), new XElement("pointsize", 1), new XElement("glyph", new XElement("visible", glyphVisible ? 1 : 0), new XElement("scale", 1)), new XElement("camera", new XElement("focumode", 0), new XElement("focus", new XAttribute("x", 0), new XAttribute("y", 0), new XAttribute("z", 0))));
		return plot;
	}

	private static XElement CreateClusterElement(int key, String label, Color color, boolean isDefault, double size, int shape)
	{
		XElement cluster = new XElement("cluster", new XElement("key", key), new XElement("label", label), new XElement("visible", 1), new XElement("default", isDefault ? 1 : 0), new XElement("color", new XAttribute("r", color.R), new XAttribute("g", color.G), new XAttribute("b", color.B), new XAttribute("a", color.A)), new XElement("size", size), new XElement("shape", shape));
		return cluster;
	}

	private static XElement CreatePointElement(int key, int clusterKey, String label, double x, double y, double z)
	{
		XElement point = new XElement("point", new XElement("key", key), new XElement("clusterkey", clusterKey), new XElement("label", label), new XElement("location", new XAttribute("x", x), new XAttribute("y", y), new XAttribute("z", z)));
		return point;
	}

//C# TO JAVA CONVERTER WARNING: Java does not allow user-defined value types. The behavior of this class will differ from the original:
//ORIGINAL LINE: struct Glyphs
	private final static class Glyphs
	{
		public static int Triangle2D = 0;
		public static int Rectangle2D = 1;
		public static int Pentagon2D = 2;
		public static int Hexagon2D = 3;
		public static int Tetrahedron3D = 4;
		public static int Cube3D = 5;
		public static int Sphere3D = 6;
		public static int Cylinder3D = 7;
	}*/
}