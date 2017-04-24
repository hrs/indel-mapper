class TabularReferenceExporter(object):
    HEADER = ["Name", "Sequence", "N20", "PAM", "Total Reads", "Cutsite",
              "Count", "Reference", "Read"]
    def __init__(self, reference_presenters):
        self.reference_presenters = reference_presenters

    def export(self):
        rows = [HEADER]
        for ref in self.reference_presenters:
            rows.extend(_rows_for(ref))
        return rows

    def _rows_for(self, ref):
        rows = []
        prefix_cells = reference_presenter.csv_row_prefix_cells()

        if reference_presenter.has_mutation_clusters():
            for cluster in reference_presenter.mutation_clusters:
                cluster_info = [cluster.cutsite_region, cluster.count()]
                for sequence in cluster.representations:
                    rows.append(prefix_cells + cluster_info +
                                [sequence.reference, sequence.read])
                if cluster.count() == 0:
                    rows.append(prefix_cells + cluster_info)
        else:
            rows.append(prefix_cells)

        return rows
