package salsa.pairwiseclustering;

public class CenterInfo {
    private int _pnum;
    private double _measure;
    private String _method;
    private int _cluster;
    private String _seqName;
    private int _seqLength;
    private String _source;
    private double _count;

    private PlotVizPoint _plotVizPoint;

    public CenterInfo(int pnum, double measure, String method, int cluster, String seqName, int seqLength) {
        _pnum = pnum;
        _measure = measure;
        _method = method;
        _cluster = cluster;
        _seqName = seqName;
        _seqLength = seqLength;
    }

    public CenterInfo(int pnum, double measure, String method, int cluster, String seqName, int seqLength,
                      String source, double count) {
        _pnum = pnum;
        _measure = measure;
        _method = method;
        _cluster = cluster;
        _seqName = seqName;
        _seqLength = seqLength;
        _source = source;
        _count = count;
    }

    public final int getPnum() {
        return _pnum;
    }

    public final void setPnum(int value) {
        _pnum = value;
    }

    public final double getMeasure() {
        return _measure;
    }

    public final void setMeasure(double value) {
        _measure = value;
    }

    public final String getMethod() {
        return _method;
    }

    public final void setMethod(String value) {
        _method = value;
    }

    public final int getCluster() {
        return _cluster;
    }

    public final void setCluster(int value) {
        _cluster = value;
    }

    public final String getSeqName() {
        return _seqName;
    }

    public final void setSeqName(String value) {
        _seqName = value;
    }

    public final int getSeqLength() {
        return _seqLength;
    }

    public final void setSeqLength(int value) {
        _seqLength = value;
    }

    public final String getSource() {
        return _source;
    }

    public final void setSource(String value) {
        _source = value;
    }

    public final double getCount() {
        return _count;
    }

    public final void setCount(double value) {
        _count = value;
    }

    public final PlotVizPoint getPlotVizPoint() {
        return _plotVizPoint;
    }

    public final void setPlotVizPoint(PlotVizPoint value) {
        _plotVizPoint = value;
    }

    public final String toString() {
        return "PointNumber=" + _pnum + "\tMeasure=" + _measure + "\tMethod=" + _method + "\tGroup=" + _cluster +
                "\tSequence=" + _seqName + "\tLength=" + _seqLength;
    }
}