import React, { useState, useEffect } from 'react'
import { ArrowLeft, Upload, FileText, Download, Calendar, User, Plus, Loader2, AlertCircle, X, File, FileSpreadsheet, Presentation, Info } from 'lucide-react'
import { useAuth } from '../hooks/useAuth'
import { supabase } from '../lib/supabase'
import { BoltBadge } from './BoltBadge'

interface ReportsPageProps {
  labId: string
  labName: string
  onBack: () => void
}

interface Report {
  id: string
  lab_id: string
  title: string
  content: string | null
  file_url: string | null
  created_by: string
  created_at: string
  creator_name?: string
}

export function ReportsPage({ labId, labName, onBack }: ReportsPageProps) {
  const { user } = useAuth()
  const [reports, setReports] = useState<Report[]>([])
  const [loading, setLoading] = useState(true)
  const [uploading, setUploading] = useState(false)
  const [error, setError] = useState('')
  const [showUploadModal, setShowUploadModal] = useState(false)
  const [showInfoModal, setShowInfoModal] = useState(false)
  const [selectedReport, setSelectedReport] = useState<Report | null>(null)
  const [uploadTitle, setUploadTitle] = useState('')
  const [uploadDescription, setUploadDescription] = useState('')
  const [selectedFile, setSelectedFile] = useState<File | null>(null)

  useEffect(() => {
    fetchReports()
  }, [labId])

  const fetchReports = async () => {
    try {
      setLoading(true)
      const { data, error } = await supabase
        .from('reports')
        .select(`
          *,
          users!reports_created_by_fkey (name)
        `)
        .eq('lab_id', labId)
        .order('created_at', { ascending: false })

      if (error) throw error

      const reportsWithCreatorNames = data?.map(report => ({
        ...report,
        creator_name: report.users?.name || 'Unknown User'
      })) || []

      setReports(reportsWithCreatorNames)
    } catch (err: any) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  const handleFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (file) {
      // Check file type
      const allowedTypes = [
        'application/pdf',
        'application/msword',
        'application/vnd.openxmlformats-officedocument.wordprocessingml.document',
        'text/plain',
        'application/vnd.ms-excel',
        'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        'application/vnd.ms-powerpoint',
        'application/vnd.openxmlformats-officedocument.presentationml.presentation'
      ]
      
      if (!allowedTypes.includes(file.type)) {
        setError('Please select a valid file type (PDF, DOC, DOCX, TXT, XLS, XLSX, PPT, PPTX)')
        return
      }

      // Check file size (50MB limit)
      if (file.size > 50 * 1024 * 1024) {
        setError('File size must be less than 50MB')
        return
      }

      setSelectedFile(file)
      setError('')
    }
  }

  const uploadReport = async () => {
    if (!selectedFile || !uploadTitle.trim() || !user) return

    setUploading(true)
    setError('')

    try {
      // Generate unique filename with user ID for better organization
      const fileExt = selectedFile.name.split('.').pop()
      const fileName = `${Date.now()}-${Math.random().toString(36).substring(2)}.${fileExt}`
      const filePath = `${labId}/${user.id}/${fileName}`

      // Upload file to storage
      const { data: uploadData, error: uploadError } = await supabase.storage
        .from('reports')
        .upload(filePath, selectedFile, {
          cacheControl: '3600',
          upsert: false
        })

      if (uploadError) {
        console.error('Upload error:', uploadError)
        throw new Error(`Upload failed: ${uploadError.message}`)
      }

      // Save report metadata to database with the storage path
      const { error: dbError } = await supabase
        .from('reports')
        .insert({
          lab_id: labId,
          title: uploadTitle.trim(),
          content: uploadDescription.trim() || null,
          file_url: filePath, // Store the storage path, not public URL
          created_by: user.id
        })

      if (dbError) throw dbError

      // Reset form
      setUploadTitle('')
      setUploadDescription('')
      setSelectedFile(null)
      setShowUploadModal(false)
      
      // Refresh reports list
      fetchReports()
    } catch (err: any) {
      console.error('Upload error:', err)
      setError(err.message)
    } finally {
      setUploading(false)
    }
  }

  const downloadReport = async (report: Report, e: React.MouseEvent) => {
    e.stopPropagation()
    if (!report.file_url) return

    try {
      // Get signed URL for download
      const { data, error } = await supabase.storage
        .from('reports')
        .createSignedUrl(report.file_url, 60) // 60 seconds expiry

      if (error) {
        console.error('Download signed URL error:', error)
        throw new Error(`Download failed: ${error.message}`)
      }

      if (data?.signedUrl) {
        // Create download link
        const link = document.createElement('a')
        link.href = data.signedUrl
        link.download = `${report.title}.${report.file_url.split('.').pop()}`
        link.target = '_blank'
        document.body.appendChild(link)
        link.click()
        document.body.removeChild(link)
      } else {
        throw new Error('No signed URL received')
      }
    } catch (err: any) {
      console.error('Download error:', err)
      setError('Failed to download file: ' + err.message)
    }
  }

  const showReportInfo = (report: Report, e: React.MouseEvent) => {
    e.stopPropagation()
    setSelectedReport(report)
    setShowInfoModal(true)
  }

  const getFileIcon = (fileUrl: string | null) => {
    if (!fileUrl) return <FileText className="w-12 h-12 text-gray-400" />
    
    const ext = fileUrl.split('.').pop()?.toLowerCase()
    switch (ext) {
      case 'pdf':
        return <FileText className="w-12 h-12 text-red-400" />
      case 'doc':
      case 'docx':
        return <File className="w-12 h-12 text-blue-400" />
      case 'xls':
      case 'xlsx':
        return <FileSpreadsheet className="w-12 h-12 text-green-400" />
      case 'ppt':
      case 'pptx':
        return <Presentation className="w-12 h-12 text-orange-400" />
      case 'txt':
        return <FileText className="w-12 h-12 text-gray-400" />
      default:
        return <FileText className="w-12 h-12 text-gray-400" />
    }
  }

  const getFileColor = (fileUrl: string | null) => {
    if (!fileUrl) return 'from-gray-800 to-gray-700'
    
    const ext = fileUrl.split('.').pop()?.toLowerCase()
    switch (ext) {
      case 'pdf':
        return 'from-red-900/30 to-red-800/30'
      case 'doc':
      case 'docx':
        return 'from-blue-900/30 to-blue-800/30'
      case 'xls':
      case 'xlsx':
        return 'from-green-900/30 to-green-800/30'
      case 'ppt':
      case 'pptx':
        return 'from-orange-900/30 to-orange-800/30'
      case 'txt':
        return 'from-gray-800 to-gray-700'
      default:
        return 'from-gray-800 to-gray-700'
    }
  }

  const getMobileFileIcon = (fileUrl: string | null) => {
    if (!fileUrl) return <FileText className="w-6 h-6 text-gray-400" />
    
    const ext = fileUrl.split('.').pop()?.toLowerCase()
    switch (ext) {
      case 'pdf':
        return <FileText className="w-6 h-6 text-red-400" />
      case 'doc':
      case 'docx':
        return <File className="w-6 h-6 text-blue-400" />
      case 'xls':
      case 'xlsx':
        return <FileSpreadsheet className="w-6 h-6 text-green-400" />
      case 'ppt':
      case 'pptx':
        return <Presentation className="w-6 h-6 text-orange-400" />
      case 'txt':
        return <FileText className="w-6 h-6 text-gray-400" />
      default:
        return <FileText className="w-6 h-6 text-gray-400" />
    }
  }

  if (loading) {
    return (
      <div className="h-screen bg-[#0F0F0F] text-white flex items-center justify-center">
        <div className="text-xl">Loading reports...</div>
      </div>
    )
  }

  return (
    <div className="h-screen bg-[#0F0F0F] text-white flex flex-col overflow-hidden">
      <div className="fixed top-[90px] right-8 z-50">
        <BoltBadge />
      </div>

      {/* Header */}
      <header className="border-b border-gray-800 bg-[#0F0F0F]/95 backdrop-blur-md flex-shrink-0">
        <div className="px-6 py-4 flex justify-between items-center">
          <div className="flex items-center gap-4">
            <button
              onClick={onBack}
              className="flex items-center gap-2 text-gray-300 hover:text-white transition-colors duration-300 group"
            >
              <ArrowLeft className="w-5 h-5 group-hover:-translate-x-1 transition-transform duration-300" />
              Back
            </button>
            <div className="h-6 w-px bg-gray-700" />
            <div className="hidden md:block">
              <h1 className="text-xl font-bold text-white">Reports</h1>
              <p className="text-sm text-gray-400">
                {labName} • {reports.length} report{reports.length !== 1 ? 's' : ''}
              </p>
            </div>
          </div>
          
          <button
            onClick={() => setShowUploadModal(true)}
            className="flex items-center gap-2 bg-blue-600 hover:bg-blue-700 text-white px-4 py-2 rounded-lg font-medium transition-colors duration-300"
          >
            <Plus className="w-4 h-4" />
            Upload Report
          </button>
        </div>
      </header>

      {/* Main Content */}
      <div className="flex-1 overflow-auto p-6">
        {/* Error Message */}
        {error && (
          <div className="mb-6 p-4 bg-red-500/20 border border-red-500/30 rounded-lg">
            <div className="flex items-center gap-2 text-red-300">
              <AlertCircle className="w-4 h-4" />
              <span className="text-sm">{error}</span>
              <button 
                onClick={() => setError('')}
                className="ml-auto text-red-400 hover:text-red-300"
              >
                <X className="w-4 h-4" />
              </button>
            </div>
          </div>
        )}

        {/* Reports Grid */}
        {reports.length === 0 ? (
          <div className="text-center py-16">
            <FileText className="w-16 h-16 mx-auto mb-4 text-gray-500" />
            <h3 className="text-xl font-semibold mb-2">No Reports Yet</h3>
            <p className="text-gray-400 mb-6">Upload your first report to get started</p>
            <button
              onClick={() => setShowUploadModal(true)}
              className="bg-blue-600 hover:bg-blue-700 text-white px-6 py-3 rounded-lg font-medium transition-colors duration-300"
            >
              Upload First Report
            </button>
          </div>
        ) : (
          <>
            {/* Desktop View - Left-aligned Compact Folders */}
            <div className="hidden md:block">
              <div className="flex flex-wrap gap-4">
                {reports.map((report) => (
                  <div
                    key={report.id}
                    className="group relative w-48"
                  >
                    {/* Dark Folder Background */}
                    <div className="relative bg-gradient-to-br from-gray-700 to-gray-800 rounded-lg p-3 shadow-lg transform transition-all duration-300 hover:scale-105 hover:shadow-xl border border-gray-600">
                      {/* Folder Tab */}
                      <div className="absolute -top-1.5 left-3 w-12 h-4 bg-gradient-to-r from-gray-600 to-gray-700 rounded-t-md border-l border-r border-t border-gray-500"></div>
                      
                      {/* Paper/Document */}
                      <div className={`bg-gradient-to-br ${getFileColor(report.file_url)} border border-gray-600 rounded-lg rounded-tl-none p-3 shadow-md h-40 flex flex-col relative`}>
                        {/* File Icon */}
                        <div className="flex justify-center mb-2">
                          {getFileIcon(report.file_url)}
                        </div>
                        
                        {/* Title on the paper */}
                        <div className="flex-1 flex items-center justify-center">
                          <h3 className="text-xs font-bold text-gray-200 text-center line-clamp-3 leading-tight px-1">
                            {report.title}
                          </h3>
                        </div>
                        
                        {/* Metadata at bottom */}
                        <div className="text-xs text-gray-400 space-y-1">
                          <div className="flex items-center justify-center gap-1">
                            <User className="w-2.5 h-2.5" />
                            <span className="truncate text-xs">{report.creator_name}</span>
                          </div>
                          <div className="flex items-center justify-center gap-1">
                            <Calendar className="w-2.5 h-2.5" />
                            <span className="text-xs">{new Date(report.created_at).toLocaleDateString()}</span>
                          </div>
                        </div>

                        {/* Hover Actions */}
                        <div className="absolute inset-0 bg-black/60 opacity-0 group-hover:opacity-100 transition-opacity duration-300 rounded-lg rounded-tl-none flex items-center justify-center">
                          <div className="flex gap-2">
                            {report.file_url && (
                              <button
                                onClick={(e) => downloadReport(report, e)}
                                className="p-2 bg-white/20 backdrop-blur-sm text-white hover:bg-white/30 rounded-full transition-all duration-300 hover:scale-110"
                                title="Download"
                              >
                                <Download className="w-4 h-4" />
                              </button>
                            )}
                            {report.content && (
                              <button
                                onClick={(e) => showReportInfo(report, e)}
                                className="p-2 bg-blue-500/20 backdrop-blur-sm text-blue-300 hover:bg-blue-500/30 rounded-full transition-all duration-300 hover:scale-110"
                                title="View Details"
                              >
                                <Info className="w-4 h-4" />
                              </button>
                            )}
                          </div>
                        </div>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>

            {/* Mobile View - Updated Layout */}
            <div className="md:hidden space-y-4">
              {reports.map((report) => (
                <div
                  key={report.id}
                  className="bg-gradient-to-br from-gray-900 to-gray-800 border border-gray-700 rounded-lg p-4"
                >
                  <div className="flex items-center justify-between">
                    {/* Left side: Icon + File Name */}
                    <div className="flex items-center gap-3 flex-1 min-w-0">
                      {getMobileFileIcon(report.file_url)}
                      <div className="flex-1 min-w-0">
                        <h3 className="text-lg font-semibold text-white truncate w-[95%]">
                          {report.title}
                        </h3>
                        <div className="flex items-center gap-3 text-xs text-gray-500 mt-1">
                          <div className="hidden md:flex items-center gap-1">
                            <User className="w-3 h-3" />
                            <span>{report.creator_name}</span>
                          </div>
                          <div className="flex items-center gap-1">
                            <Calendar className="w-3 h-3" />
                            <span>{new Date(report.created_at).toLocaleDateString()}</span>
                          </div>
                        </div>
                      </div>
                    </div>
                    
                    {/* Right side: Description + Info Icon + Download Button */}
                    <div className="flex items-center gap-2 flex-shrink-0">
                      {/* Description and Info Icon - only show if description exists */}
                      {report.content && (
                        <button
                          onClick={(e) => showReportInfo(report, e)}
                          className="p-2 bg-blue-500/20 text-blue-300 hover:bg-blue-500/30 rounded-full transition-all duration-300"
                          title="View Details"
                        >
                          <Info className="w-4 h-4" />
                        </button>
                      )}
                      
                      {/* Download Button */}
                      {report.file_url && (
                        <button
                          onClick={(e) => downloadReport(report, e)}
                          className="flex items-center gap-1 px-3 py-2 bg-blue-600 hover:bg-blue-700 text-white text-sm rounded-lg transition-colors duration-300"
                        >
                          <Download className="w-4 h-4" />
                          Download
                        </button>
                      )}
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </>
        )}
      </div>

      {/* Upload Modal */}
      {showUploadModal && (
        <div className="fixed inset-0 z-50 overflow-y-auto">
          <div className="fixed inset-0 bg-black/60 backdrop-blur-sm" onClick={() => setShowUploadModal(false)} />
          <div className="flex min-h-full items-center justify-center p-4">
            <div className="relative w-full max-w-md transform overflow-hidden rounded-2xl bg-gradient-to-br from-gray-900 to-gray-800 border border-gray-700 p-8 shadow-2xl transition-all">
              <button
                onClick={() => {
                  setShowUploadModal(false)
                  setUploadTitle('')
                  setUploadDescription('')
                  setSelectedFile(null)
                  setError('')
                }}
                className="absolute right-4 top-4 text-gray-400 hover:text-white transition-colors duration-300"
              >
                <X className="w-6 h-6" />
              </button>

              <div className="mb-6">
                <div className="bg-gradient-to-r from-blue-500 to-purple-600 w-14 h-14 rounded-xl flex items-center justify-center mb-4">
                  <Upload className="w-7 h-7 text-white" />
                </div>
                <h2 className="text-2xl font-bold text-white mb-2">Upload Report</h2>
                <p className="text-gray-400">Add a new report to your lab</p>
              </div>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-gray-300 mb-2">
                    Report Title *
                  </label>
                  <input
                    type="text"
                    value={uploadTitle}
                    onChange={(e) => setUploadTitle(e.target.value)}
                    className="w-full px-4 py-3 bg-gray-800 border border-gray-600 rounded-lg text-white placeholder-gray-400 focus:outline-none focus:ring-2 focus:ring-blue-500"
                    placeholder="Enter report title"
                    disabled={uploading}
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium text-gray-300 mb-2">
                    Description
                  </label>
                  <textarea
                    value={uploadDescription}
                    onChange={(e) => setUploadDescription(e.target.value)}
                    rows={3}
                    className="w-full px-4 py-3 bg-gray-800 border border-gray-600 rounded-lg text-white placeholder-gray-400 focus:outline-none focus:ring-2 focus:ring-blue-500 resize-none"
                    placeholder="Optional description"
                    disabled={uploading}
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium text-gray-300 mb-2">
                    File *
                  </label>
                  <input
                    type="file"
                    onChange={handleFileSelect}
                    accept=".pdf,.doc,.docx,.txt,.xls,.xlsx,.ppt,.pptx"
                    className="w-full px-4 py-3 bg-gray-800 border border-gray-600 rounded-lg text-white file:mr-4 file:py-2 file:px-4 file:rounded file:border-0 file:text-sm file:font-medium file:bg-blue-600 file:text-white hover:file:bg-blue-700"
                    disabled={uploading}
                  />
                  {selectedFile && (
                    <div className="mt-2 text-sm text-gray-400">
                      Selected: {selectedFile.name} ({(selectedFile.size / 1024 / 1024).toFixed(2)} MB)
                    </div>
                  )}
                </div>

                <div className="flex gap-4 pt-4">
                  <button
                    onClick={() => {
                      setShowUploadModal(false)
                      setUploadTitle('')
                      setUploadDescription('')
                      setSelectedFile(null)
                      setError('')
                    }}
                    disabled={uploading}
                    className="flex-1 px-6 py-3 border border-gray-600 text-gray-300 rounded-lg font-medium hover:border-gray-500 hover:text-white transition-all duration-300 disabled:opacity-50 disabled:cursor-not-allowed"
                  >
                    Cancel
                  </button>
                  <button
                    onClick={uploadReport}
                    disabled={uploading || !uploadTitle.trim() || !selectedFile}
                    className="flex-1 bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 disabled:from-gray-600 disabled:to-gray-700 text-white px-6 py-3 rounded-lg font-medium transition-all duration-300 disabled:cursor-not-allowed flex items-center justify-center gap-2"
                  >
                    {uploading ? (
                      <>
                        <Loader2 className="w-4 h-4 animate-spin" />
                        Uploading...
                      </>
                    ) : (
                      <>
                        <Upload className="w-4 h-4" />
                        Upload
                      </>
                    )}
                  </button>
                </div>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Info Modal */}
      {showInfoModal && selectedReport && (
        <div className="fixed inset-0 z-50 overflow-y-auto">
          <div className="fixed inset-0 bg-black/60 backdrop-blur-sm" onClick={() => setShowInfoModal(false)} />
          <div className="flex min-h-full items-center justify-center p-4">
            <div className="relative w-full max-w-md transform overflow-hidden rounded-2xl bg-gradient-to-br from-gray-900 to-gray-800 border border-gray-700 p-8 shadow-2xl transition-all">
              <button
                onClick={() => setShowInfoModal(false)}
                className="absolute right-4 top-4 text-gray-400 hover:text-white transition-colors duration-300"
              >
                <X className="w-6 h-6" />
              </button>

              <div className="mb-6">
                <div className="bg-gradient-to-r from-blue-500 to-cyan-600 w-14 h-14 rounded-xl flex items-center justify-center mb-4">
                  <Info className="w-7 h-7 text-white" />
                </div>
                <h2 className="text-2xl font-bold text-white mb-2">Report Details</h2>
              </div>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-gray-400 mb-1">Title</label>
                  <p className="text-white font-medium">{selectedReport.title}</p>
                </div>

                {selectedReport.content && (
                  <div>
                    <label className="block text-sm font-medium text-gray-400 mb-1">Description</label>
                    <p className="text-gray-300 leading-relaxed">{selectedReport.content}</p>
                  </div>
                )}

                <div>
                  <label className="block text-sm font-medium text-gray-400 mb-1">Created By</label>
                  <p className="text-gray-300">{selectedReport.creator_name}</p>
                </div>

                <div>
                  <label className="block text-sm font-medium text-gray-400 mb-1">Created On</label>
                  <p className="text-gray-300">{new Date(selectedReport.created_at).toLocaleDateString('en-US', {
                    year: 'numeric',
                    month: 'long',
                    day: 'numeric',
                    hour: '2-digit',
                    minute: '2-digit'
                  })}</p>
                </div>

                {selectedReport.file_url && (
                  <div className="pt-4">
                    <button
                      onClick={(e) => {
                        downloadReport(selectedReport, e)
                        setShowInfoModal(false)
                      }}
                      className="w-full bg-gradient-to-r from-blue-600 to-cyan-600 hover:from-blue-700 hover:to-cyan-700 text-white px-6 py-3 rounded-lg font-medium transition-all duration-300 flex items-center justify-center gap-2"
                    >
                      <Download className="w-4 h-4" />
                      Download Report
                    </button>
                  </div>
                )}
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}