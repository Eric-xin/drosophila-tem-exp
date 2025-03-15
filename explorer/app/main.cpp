#include <QApplication>
#include <QMainWindow>
#include <QTableView>
#include <QAbstractTableModel>
#include <QFileDialog>
#include <QComboBox>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QDialog>
#include <QGridLayout>
#include <QTextEdit>
#include <QProgressDialog>
#include <QMessageBox>
#include <QGroupBox>
#include <QScrollArea>
#include <QClipboard>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QUrlQuery>
#include <QEventLoop>
#include <QThread>
#include <QFile>
#include <QTextStream>
#include <QSet>
#include <QVariantMap>
#include <algorithm>

// ---------------------------
// TableModel: holds TSV data
// ---------------------------
class TableModel : public QAbstractTableModel {
    Q_OBJECT
public:
    TableModel(QObject *parent = nullptr) : QAbstractTableModel(parent) {}
    
    int rowCount(const QModelIndex &parent = QModelIndex()) const override {
        Q_UNUSED(parent);
        return m_data.size();
    }
    int columnCount(const QModelIndex &parent = QModelIndex()) const override {
        Q_UNUSED(parent);
        return m_headers.size();
    }
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override {
        if (!index.isValid())
            return QVariant();
        const QVariantMap &row = m_data.at(index.row());
        QString key = m_headers.at(index.column());
        if (role == Qt::DisplayRole) {
            QVariant value = row.value(key);
            if (value.isNull())
                return "N/A";
            if (key == "start" || key == "end")
                return QString::number(value.toInt());
            if (value.type() == QVariant::Double)
                return QString::number(value.toDouble(), 'f', 4);
            return value.toString();
        }
        if (role == Qt::TextAlignmentRole) {
            QVariant value = row.value(key); // redeclare value here
            if (value.canConvert<double>() && !value.isNull())
                return QVariant(int(Qt::AlignRight | Qt::AlignVCenter));
            else
                return QVariant(int(Qt::AlignLeft | Qt::AlignVCenter));
        }
        return QVariant();
    }
    QVariant headerData(int section, Qt::Orientation orientation, int role=Qt::DisplayRole) const override {
        if (role != Qt::DisplayRole)
            return QVariant();
        if (orientation == Qt::Horizontal) {
            if (section < m_headers.size())
                return m_headers.at(section);
            return QVariant();
        } else {
            return QString::number(section);
        }
    }
    void updateData(const QList<QVariantMap> &data, const QStringList &headers) {
        beginResetModel();
        m_data = data;
        m_headers = headers;
        endResetModel();
    }
    QList<QVariantMap> getData() const { return m_data; }
private:
    QList<QVariantMap> m_data;
    QStringList m_headers;
};

// --------------------------------------
// Utility: Load a TSV file into memory
// --------------------------------------
bool loadTSV(const QString &fileName, QList<QVariantMap> &data, QStringList &headers) {
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return false;
    QTextStream in(&file);
    bool firstLine = true;
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList parts = line.split("\t");
        if (firstLine) {
            headers = parts;
            firstLine = false;
        } else {
            QVariantMap row;
            for (int i = 0; i < headers.size() && i < parts.size(); i++) {
                bool ok;
                double d = parts[i].toDouble(&ok);
                if (ok)
                    row[headers.at(i)] = d;
                else
                    row[headers.at(i)] = parts.at(i);
            }
            data.append(row);
        }
    }
    file.close();
    return true;
}

// -----------------------------------------------------
// Functions to fetch external data using FlyBase API
// (Synchronous calls using QNetworkAccessManager)
// -----------------------------------------------------
QString fetchSequenceId(const QString &geneId) {
    QUrl url(QString("https://api.flybase.org/api/v1.0/sequence/id/%1").arg(geneId));
    QNetworkRequest request(url);
    request.setRawHeader("accept", "application/json");
    QNetworkAccessManager manager;
    QEventLoop loop;
    QNetworkReply *reply = manager.get(request);
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();
    if (reply->error() == QNetworkReply::NoError) {
        QByteArray response_data = reply->readAll();
        QJsonDocument doc = QJsonDocument::fromJson(response_data);
        QJsonObject obj = doc.object();
        if (obj.contains("resultset")) {
            QJsonObject rs = obj["resultset"].toObject();
            if (rs.contains("result")) {
                QJsonArray arr = rs["result"].toArray();
                if (!arr.isEmpty()) {
                    QJsonObject first = arr.first().toObject();
                    return first["sequence"].toString();
                }
            }
        }
        return "No sequence found for the given gene ID.";
    } else {
        return QString("Error: %1").arg(reply->errorString());
    }
}

QString fetchSequenceLocation(const QString &species, const QString &location, const QString &strand="minus", int padding=100) {
    QString loc = location;
    if (loc.contains(".."))
        loc.replace("..", "-");
    QUrl url(QString("https://api.flybase.org/api/v1.0/sequence/region/%1/%2").arg(species, loc));
    QUrlQuery query;
    query.addQueryItem("strand", strand);
    query.addQueryItem("padding", QString::number(padding));
    url.setQuery(query);
    QNetworkRequest request(url);
    request.setRawHeader("accept", "application/json");
    QNetworkAccessManager manager;
    QEventLoop loop;
    QNetworkReply *reply = manager.get(request);
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();
    if (reply->error() == QNetworkReply::NoError) {
        QByteArray response_data = reply->readAll();
        QJsonDocument doc = QJsonDocument::fromJson(response_data);
        QJsonObject obj = doc.object();
        if (obj.contains("resultset")) {
            QJsonObject rs = obj["resultset"].toObject();
            if (rs.contains("result")) {
                QJsonArray arr = rs["result"].toArray();
                if (!arr.isEmpty()) {
                    QJsonObject first = arr.first().toObject();
                    return first["sequence"].toString();
                }
            }
        }
        return "No sequence returned for this region.";
    } else {
        return QString("Error: %1").arg(reply->errorString());
    }
}

QString fetchIdAbt(const QString &geneId) {
    QUrl url(QString("https://api.flybase.org/api/v1.0/gene/summaries/auto/%1").arg(geneId));
    QNetworkRequest request(url);
    request.setRawHeader("accept", "application/json");
    QNetworkAccessManager manager;
    QEventLoop loop;
    QNetworkReply *reply = manager.get(request);
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();
    if (reply->error() == QNetworkReply::NoError) {
        QByteArray response_data = reply->readAll();
        QJsonDocument doc = QJsonDocument::fromJson(response_data);
        QJsonObject obj = doc.object();
        if (obj.contains("resultset")) {
            QJsonObject rs = obj["resultset"].toObject();
            if (rs.contains("result")) {
                QJsonArray arr = rs["result"].toArray();
                if (!arr.isEmpty()) {
                    QJsonObject first = arr.first().toObject();
                    return first["summary"].toString();
                }
            }
        }
        return "No summary found for the given gene ID.";
    } else {
        return QString("Error: %1").arg(reply->errorString());
    }
}

// ---------------------------------------
// Utility: Format sequence with HTML tags
// ---------------------------------------
QString formatSequence(const QString &seq) {
    QString html = "<div style=\"white-space: normal;\">";
    for (const QChar &ch : seq) {
        QString up = ch.toUpper();
        if (up == "A")
            html += "<span style=\"color: red; font-weight: bold;\">" + QString(ch) + "</span>";
        else if (up == "T")
            html += "<span style=\"color: blue; font-weight: bold;\">" + QString(ch) + "</span>";
        else if (up == "C")
            html += "<span style=\"color: green; font-weight: bold;\">" + QString(ch) + "</span>";
        else if (up == "G")
            html += "<span style=\"color: orange; font-weight: bold;\">" + QString(ch) + "</span>";
        else
            html += ch;
    }
    html += "</div>";
    return html;
}

// -----------------------------------------
// DetailDialog: shows details for a gene
// -----------------------------------------
class DetailDialog : public QDialog {
    Q_OBJECT
public:
    DetailDialog(const QVariantMap &rowData, QWidget *parent = nullptr) : QDialog(parent), m_details(rowData) {
        setWindowTitle("Gene Detail");
        resize(700, 500);
        QVBoxLayout *mainLayout = new QVBoxLayout(this);
        
        QPushButton *copyButton = new QPushButton("Copy as JSON", this);
        connect(copyButton, &QPushButton::clicked, this, &DetailDialog::copyAsJson);
        mainLayout->addWidget(copyButton);
        
        QScrollArea *scroll = new QScrollArea(this);
        scroll->setWidgetResizable(true);
        mainLayout->addWidget(scroll);
        QWidget *content = new QWidget();
        scroll->setWidget(content);
        QVBoxLayout *contentLayout = new QVBoxLayout(content);
        
        // Basic fields in a grid (excluding fbgn, arm, start, end)
        QGridLayout *grid = new QGridLayout();
        int idx = 0;
        QSet<QString> excludeKeys = {"fbgn", "arm", "start", "end"};
        for (auto it = m_details.begin(); it != m_details.end(); ++it) {
            if (excludeKeys.contains(it.key().toLower()))
                continue;
            QGroupBox *box = createFieldBox(it.key(), it.value().toString());
            grid->addWidget(box, idx / 2, idx % 2);
            idx++;
        }
        contentLayout->addLayout(grid);
        
        // External details: Gene Summary and Gene Sequence
        QString fbgn = m_details.value("fbgn").toString().trimmed();
        QString arm = m_details.value("arm").toString().trimmed();
        QString geneSummary, geneSequence;
        if (!fbgn.isEmpty() && fbgn.startsWith("FBgn")) {
            geneSummary = fetchIdAbt(fbgn);
            QGroupBox *summaryBox = createFieldBox("Gene Summary", geneSummary);
            contentLayout->addWidget(summaryBox);
            geneSequence = fetchSequenceId(fbgn);
            QGroupBox *seqBox = createSequenceBox("Gene Sequence", geneSequence);
            contentLayout->addWidget(seqBox);
        } else if (!arm.isEmpty() && m_details.contains("start") && m_details.contains("end")) {
            QString location = QString("%1:%2..%3").arg(arm)
                                .arg(m_details.value("start").toString())
                                .arg(m_details.value("end").toString());
            geneSequence = fetchSequenceLocation("dmel", location);
            QGroupBox *seqBox = createSequenceBox(QString("Gene Sequence (from %1)").arg(location), geneSequence);
            contentLayout->addWidget(seqBox);
        }
        // Store external details
        m_details["Gene Summary"] = geneSummary;
        m_details["Gene Sequence"] = geneSequence;
    }
private slots:
    void copyAsJson() {
        QJsonObject jsonObj;
        for (auto it = m_details.begin(); it != m_details.end(); ++it)
            jsonObj.insert(it.key(), QJsonValue::fromVariant(it.value()));
        QJsonDocument doc(jsonObj);
        QString jsonString = doc.toJson(QJsonDocument::Indented);
        QApplication::clipboard()->setText(jsonString);
        QMessageBox::information(this, "Copied", "Data copied as JSON to clipboard.");
    }
private:
    QGroupBox* createFieldBox(const QString &title, const QString &content) {
        QGroupBox *box = new QGroupBox(title);
        QVBoxLayout *layout = new QVBoxLayout(box);
        QTextEdit *textEdit = new QTextEdit();
        textEdit->setPlainText(content);
        textEdit->setReadOnly(true);
        textEdit->setFrameStyle(QFrame::NoFrame);
        textEdit->setLineWrapMode(QTextEdit::WidgetWidth);
        layout->addWidget(textEdit);
        return box;
    }
    QGroupBox* createSequenceBox(const QString &title, const QString &sequence) {
        QGroupBox *box = new QGroupBox(title);
        QVBoxLayout *layout = new QVBoxLayout(box);
        QTextEdit *textEdit = new QTextEdit();
        textEdit->setReadOnly(true);
        textEdit->setFrameStyle(QFrame::NoFrame);
        textEdit->setLineWrapMode(QTextEdit::WidgetWidth);
        if (!sequence.isEmpty())
            textEdit->setHtml(formatSequence(sequence));
        else
            textEdit->setPlainText("No sequence available.");
        layout->addWidget(textEdit);
        return box;
    }
    QVariantMap m_details;
};

// -----------------------------------------------------
// DataExplorer: Main window with table and controls
// -----------------------------------------------------
class DataExplorer : public QMainWindow {
    Q_OBJECT
public:
    DataExplorer(const QString &preloadFile = QString(), QWidget *parent = nullptr) : QMainWindow(parent) {
        setWindowTitle("Data Explorer");
        resize(1000, 600);
        m_model = new TableModel(this);
        
        QWidget *mainWidget = new QWidget(this);
        setCentralWidget(mainWidget);
        QVBoxLayout *mainLayout = new QVBoxLayout(mainWidget);
        
        // Top control panel
        QHBoxLayout *controlLayout = new QHBoxLayout();
        QPushButton *loadButton = new QPushButton("Load TSV Data", this);
        connect(loadButton, &QPushButton::clicked, this, &DataExplorer::openFileDialog);
        controlLayout->addWidget(loadButton);
        
        QLabel *sortLabel = new QLabel("Sort Mode:", this);
        controlLayout->addWidget(sortLabel);
        
        m_sortModeCombo = new QComboBox(this);
        m_sortModeCombo->addItems({"Sort by Location", "Sort by Temperature Reading", "Sort by Temperature Ratio"});
        connect(m_sortModeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &DataExplorer::updateTemperatureControls);
        controlLayout->addWidget(m_sortModeCombo);
        
        m_tempSortLabel = new QLabel("Temperature Reading:", this);
        m_tempSortCombo = new QComboBox(this);
        m_tempSortCombo->addItems({"temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"});
        controlLayout->addWidget(m_tempSortLabel);
        controlLayout->addWidget(m_tempSortCombo);
        
        m_ratioNumLabel = new QLabel("Temperature Numerator:", this);
        m_ratioNumCombo = new QComboBox(this);
        m_ratioNumCombo->addItems({"temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"});
        m_ratioDenomLabel = new QLabel("Temperature Denom:", this);
        m_ratioDenomCombo = new QComboBox(this);
        m_ratioDenomCombo->addItems({"temp13_avg", "temp18_avg", "temp23_avg", "temp29_avg"});
        controlLayout->addWidget(m_ratioNumLabel);
        controlLayout->addWidget(m_ratioNumCombo);
        controlLayout->addWidget(m_ratioDenomLabel);
        controlLayout->addWidget(m_ratioDenomCombo);
        
        QPushButton *sortButton = new QPushButton("Sort Data", this);
        connect(sortButton, &QPushButton::clicked, this, &DataExplorer::sortData);
        controlLayout->addWidget(sortButton);
        mainLayout->addLayout(controlLayout);
        
        updateTemperatureControls();
        
        m_tableView = new QTableView(this);
        m_tableView->setModel(m_model);
        m_tableView->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        connect(m_tableView, &QTableView::doubleClicked, this, &DataExplorer::tableDoubleClick);
        mainLayout->addWidget(m_tableView);
        
        if (!preloadFile.isEmpty())
            loadFile(preloadFile);
    }
private slots:
    void openFileDialog() {
        QString filename = QFileDialog::getOpenFileName(this, "Open TSV File", "", "TSV Files (*.tsv *.txt);;All Files (*)");
        if (!filename.isEmpty())
            loadFile(filename);
    }
    void loadFile(const QString &filename) {
        QList<QVariantMap> data;
        QStringList headers;
        if (!loadTSV(filename, data, headers)) {
            QMessageBox::critical(this, "Error", "Failed to load file.");
            return;
        }
        m_allData = data;
        m_headers = headers;
        m_model->updateData(data, headers);
    }
    void sortData() {
        if (m_allData.isEmpty())
            return;
        QString mode = m_sortModeCombo->currentText();
        QList<QVariantMap> sortedData = m_allData; // copy data
        if (mode == "Sort by Location") {
            std::sort(sortedData.begin(), sortedData.end(), [](const QVariantMap &a, const QVariantMap &b) {
                QString armA = a.value("arm").toString();
                QString armB = b.value("arm").toString();
                if (armA != armB)
                    return armA < armB;
                int startA = a.value("start").toInt();
                int startB = b.value("start").toInt();
                if (startA != startB)
                    return startA < startB;
                return a.value("end").toInt() < b.value("end").toInt();
            });
        } else if (mode == "Sort by Temperature Reading") {
            QString sortCol = m_tempSortCombo->currentText();
            std::sort(sortedData.begin(), sortedData.end(), [sortCol](const QVariantMap &a, const QVariantMap &b) {
                return a.value(sortCol).toDouble() > b.value(sortCol).toDouble();
            });
        } else if (mode == "Sort by Temperature Ratio") {
            QString numCol = m_ratioNumCombo->currentText();
            QString denomCol = m_ratioDenomCombo->currentText();
            if (numCol == denomCol) {
                QMessageBox::critical(this, "Error", "Please select two different temperature columns for the ratio.");
                return;
            }
            std::sort(sortedData.begin(), sortedData.end(), [numCol, denomCol](const QVariantMap &a, const QVariantMap &b) {
                double aNum = a.value(numCol).toDouble();
                double aDenom = a.value(denomCol).toDouble();
                double bNum = b.value(numCol).toDouble();
                double bDenom = b.value(denomCol).toDouble();
                double ratioA = (aDenom == 0 ? 0 : aNum / aDenom);
                double ratioB = (bDenom == 0 ? 0 : bNum / bDenom);
                return ratioA > ratioB;
            });
        }
        m_allData = sortedData;
        m_model->updateData(sortedData, m_headers);
    }
    void tableDoubleClick(const QModelIndex &index) {
        if (!index.isValid())
            return;
        int row = index.row();
        if (row < 0 || row >= m_allData.size())
            return;
        QVariantMap rowData = m_allData.at(row);
        DetailDialog dialog(rowData, this);
        dialog.exec();
    }
    void updateTemperatureControls() {
        QString mode = m_sortModeCombo->currentText();
        bool showTemp = (mode == "Sort by Temperature Reading");
        bool showRatio = (mode == "Sort by Temperature Ratio");
        m_tempSortLabel->setVisible(showTemp);
        m_tempSortCombo->setVisible(showTemp);
        m_ratioNumLabel->setVisible(showRatio);
        m_ratioNumCombo->setVisible(showRatio);
        m_ratioDenomLabel->setVisible(showRatio);
        m_ratioDenomCombo->setVisible(showRatio);
    }
private:
    TableModel *m_model;
    QTableView *m_tableView;
    QComboBox *m_sortModeCombo;
    QLabel *m_tempSortLabel;
    QComboBox *m_tempSortCombo;
    QLabel *m_ratioNumLabel;
    QComboBox *m_ratioNumCombo;
    QLabel *m_ratioDenomLabel;
    QComboBox *m_ratioDenomCombo;
    QList<QVariantMap> m_allData;
    QStringList m_headers;
};

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    QString preloadFile = ":/data/data.tsv";
    if (argc > 1)
        preloadFile = argv[1];
    DataExplorer explorer(preloadFile);
    explorer.show();
    return app.exec();
}

#include "main.moc"